%% Plot scatters for pairwise correlations
function showScatterPlots(xVal,medPairCorrR,medCorrEnvelopeR,medCorrInfraSlowR,...
    xLabel,yLabel,xLim,yLim,textLocX,textLocY1,textLocY2,bandLabels)
figure;
plotIdx = 1;
for iType = 1:3
    switch iType
        case 1
            plotVal   = medPairCorrR;
            plotLabel = 'Timeseries';
        case 2
            plotVal   = medCorrEnvelopeR;
            plotLabel = 'Envelope';
        case 3
            plotVal   = medCorrInfraSlowR;
            plotLabel = 'Infraslow';
    end

    for iBand = 1:5
        subplot(3,5,plotIdx);
        if strcmp(xLabel,'Distance')
            showExpFit(xVal,plotVal(:,iBand),textLocX,textLocY1,textLocY2)
        else
            showLinearFit(xVal,plotVal(:,iBand),textLocX,textLocY1,textLocY2)
             % showLinearFit(plotVal(:,iBand),xVal,textLocX,textLocY1,textLocY2)
        end
        xlim(xLim); ylim(yLim); box off;
        xlabel(xLabel); ylabel(yLabel);
        title([plotLabel '-' bandLabels{iBand}]);
        plotIdx = plotIdx+1;
    end
end
end

%% Function to fit a line
function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
coeff = polyfit(xVal,yVal,1);
xFit  = linspace(min(xVal),max(xVal),1000);
yFit  = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
plot(xFit,yFit,'-k','LineWidth',1);
text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end

%% Fit exponential curve
function showExpFit(xVal,yVal,textLocX,textLocY1,textLocY2)
plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;

% Fit exponential function
options  = optimoptions('lsqcurvefit', 'Display', 'off','Algorithm','levenberg-marquardt');
modelfun = @(b,x) b(1) * exp(-b(2).*x);
x0       = double([1 mean(yVal,'omitnan')]); % Set initial values to mean of x for better estimation of model parameters
beta0    = lsqcurvefit(modelfun,x0,xVal,double(yVal),[],[],options); % Optimize initial values 

mdl = fitnlm(xVal,yVal, modelfun, beta0);
X   = linspace(min(xVal),max(xVal),1000);

coefficients = mdl.Coefficients{:, 'Estimate'};
yFitted      = coefficients(1) * exp(-coefficients(2).*X);

plot(X,yFitted, '-k', 'LineWidth',1);
text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end