monkeys   = {'CharlieSheen'; 'Whiskey'};
commonDir = 'C:\Users\kem294\Documents\Data';

cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\DualProbe\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\DualProbe\neuroshare']));
addpath(genpath([commonDir '\Codes\DualProbe']));
addpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12\spectral_analysis\continuous\dupes']));

hemisphere     = 'Left';

% Combine all monkey data
for iM = 1:2
    allMonkeyVars(iM) =  load(['D:\Data\' monkeys{iM} '_SqM\Left Hemisphere\DualProbeVars.mat']);
end

connVals = [allMonkeyVars(1).connValsR; allMonkeyVars(2).connValsR];
distVals = [allMonkeyVars(1).distValsR; allMonkeyVars(2).distValsR];

% Pairwise correlations
medPairCorr      = [allMonkeyVars(1).medPairCorrR; allMonkeyVars(2).medPairCorrR];
medEnvelopeCorr  = [allMonkeyVars(1).medCorrEnvelopeR; allMonkeyVars(2).medCorrEnvelopeR];
medInfraSlowCorr = [allMonkeyVars(1).medCorrInfraSlowR; allMonkeyVars(2).medCorrInfraSlowR];
%%
showScatterPlots(connVals,medPairCorr,medEnvelopeCorr,medInfraSlowCorr,...
    'Functional Connectivity','Pairwise correlations',[-0.6 1],[-0.6 1],...
    0.8,-0.5,-0.4,bandLabels)

%% Pairwise correlations vs distance
showScatterPlots(distVals,medPairCorr,medEnvelopeCorr,medInfraSlowCorr,'Distance',...
    'Pairwise correlations',[0 20],[-0.6 1],12,-0.5,-0.4,bandLabels);

%% Distributions of pairwise correlations and distance
figure; subplot(121); histogram(connVals,-0.6:0.2:1); box off;xlabel('Functional connectivity');ylim([0 50]);xticks(-0.6:0.2:1);
subplot(122); histogram(distVals,0:2:16); xlabel('Distance (mm)'); box off;ylim([0 50]); xticks(0:2:20);

%% Partial correlation calculations
for iType = 1:3
    clear var
    switch iType
        case 1
            var = medPairCorr;
        case 2
            var = medEnvelopeCorr;
        case 3
            var = medInfraSlowCorr; 
    end
    for iBand = 1:5
        [rhoVal(iType,iBand),pVal(iType,iBand)] = partialcorr(connVals,var(:,iBand),distVals,"Type","Pearson");
    end
end

%% Bar plot of correlations vs timescales
figure;
bar(rhoVal); legend(bandLabels,'Location','southeast'); xticklabels({'Time series';'Power';'Infraslow'});

%% Bar plot of correlations for power
figure; bar(rhoVal(2,[2 4]));xticklabels({'Alpha','Gamma'});box off;
%% Plot pairwise correlations, FC and distance
figure; scatter3(distVals,medEnvelopeCorr(:,4),connVals,'filled');
xlabel('Distance (mm)'); ylabel('Power'); zlabel('FC');
xlim([0 20]); ylim([-0.6 1]); zlim([-0.6 1]);

mdl = fitlm([medEnvelopeCorr(:,4) distVals],connVals);
mdl2 = fitlm([medEnvelopeCorr(:,4) distVals],connVals,'y ~1 + x1 + x2^2');
mdl3 = fitlm([medEnvelopeCorr(:,3) medEnvelopeCorr(:,4) distVals],connVals);
mdl4 = fitlm([medEnvelopeCorr(:,2) medEnvelopeCorr(:,3) distVals],connVals);

envNorm  = (medEnvelopeCorr(:,4) - mean(medEnvelopeCorr(:,4)))./std(medEnvelopeCorr(:,4));
distNorm = (distVals - mean(distVals))./std(distVals);
mdlNorm  = fitlm([envNorm distNorm],connVals);

% Ranking the co-efficients

b0 = table2array(mdlNorm.Coefficients(1,1));
b1 = table2array(mdlNorm.Coefficients(2,1));
b2 = table2array(mdlNorm.Coefficients(3,1));

ci0 = [b0-1.96*table2array(mdlNorm.Coefficients(1,2)); b0+1.96*table2array(mdlNorm.Coefficients(1,2))];
ci1 = [b0-1.96*table2array(mdlNorm.Coefficients(2,2)); b1+1.96*table2array(mdlNorm.Coefficients(2,2))];
ci2 = [b0-1.96*table2array(mdlNorm.Coefficients(3,2)); b2+1.96*table2array(mdlNorm.Coefficients(3,2))];


figure; plot([b0 b1 b2],[1 2 3],'k.','MarkerSize', 25);
ylim([0 4]); xlim([-0.5 1]); box off; yticks(1:3); yticklabels(mdlNorm.CoefficientNames);
xline(0); hold on; xticks(-0.5:0.2:1);
figure; errorbar([b0 b1 b2],[1 2 3],[ci0(1) ci1(1) ci2(1)],[ci0(2) ci1(2) ci2(2)],'.','horizontal',"MarkerSize",10, "MarkerEdgeColor","k","MarkerFaceColor",'k');

%% Checking fit
yVal = medEnvelopeCorr(:,4); xVal = connVals;
figure; plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
coeff = polyfit(xVal,yVal,2);
xFit  = linspace(min(xVal),max(xVal),1000);
yFit  = polyval(coeff,xFit); mdl = fitlm(xVal,yVal,'quadratic');
plot(xFit,yFit,'-k','LineWidth',1);
text(-0.2, 0.7,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(-0.2, 0.65,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
