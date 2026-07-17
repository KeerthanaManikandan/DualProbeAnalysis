%% combinedAnalysisEphysDualProbeRS_vFinal
% Combined analysis for both animals
% Keerthana Manikandan

commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\DualProbe\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\DualProbe\neuroshare']));
addpath(genpath([commonDir '\Codes\DualProbe']));
addpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12\spectral_analysis\continuous\dupes']));

%% Combine all monkey data and initialize variables
clear allMonkeyVars
monkeys     = {'CharlieSheen'; 'Whiskey'};
bandLabels  = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking','0-30 Hz'};
timeLabels  = {'Time series','Power','Infraslow'};
layerLabels = {'S','M','D'};
areaLabels  = {'S-S','M-M','S-M'};
hemisphere  = 'Left';

for iM = 1:2
    allMonkeyVars(iM) =  load(['D:\Data\' monkeys{iM} '_SqM\Left Hemisphere\DualProbeVarsBinned.mat']);
end

connVals = [allMonkeyVars(1).connValsR; allMonkeyVars(2).connValsR];
distVals = [allMonkeyVars(1).distValsR; allMonkeyVars(2).distValsR];

% Pairwise correlations
medPairCorr      = [allMonkeyVars(1).medPairCorrR; allMonkeyVars(2).medPairCorrR];
medEnvelopeCorr  = [allMonkeyVars(1).medCorrEnvelopeR; allMonkeyVars(2).medCorrEnvelopeR];
medInfraSlowCorr = [allMonkeyVars(1).medCorrInfraSlowR; allMonkeyVars(2).medCorrInfraSlowR];


% Get cortical areas of the pairs
pairClass = [allMonkeyVars(1).pairClass; allMonkeyVars(2).pairClass];
smLoc = sum(pairClass=='SM',2)==2;
ssLoc = sum(pairClass=='SS',2)==2;
mmLoc = sum(pairClass=='MM',2)==2;


connValsZ  = zscore(connVals,[],1);
medCorrZ   = zscore(medPairCorr,[],1);
envCorrZ   = zscore(medEnvelopeCorr,[],1); 
infraCorrZ = zscore(medInfraSlowCorr,[],1) ;


% Rescaling
maxVal = 1; minVal = 0;
connValsZ  =  ((maxVal-minVal).*(connValsZ-min(connValsZ))./(max(connValsZ)-min(connValsZ)))+minVal;
medCorrZ   =  ((maxVal-minVal).*(medCorrZ-min(medCorrZ))./(max(medCorrZ)-min(medCorrZ)))+minVal;
envCorrZ   =  ((maxVal-minVal).*(envCorrZ-min(envCorrZ))./(max(envCorrZ)-min(envCorrZ)))+minVal;
infraCorrZ =  ((maxVal-minVal).*(infraCorrZ-min(infraCorrZ))./(max(infraCorrZ)-min(infraCorrZ)))+minVal;

showScatterPlots(connValsZ,medCorrZ,envCorrZ,infraCorrZ,...
'Functional Connectivity','Pairwise correlations',[-1 1],[-1 1],...
-0.5,0.9,0.8,bandLabels);

% Rescale between -1 to 1

% Normalize FC and correlations
devPairCorr = abs(medCorrZ-connValsZ);
devEnvelope = abs(envCorrZ-connValsZ);
devInfraSlow = abs(infraCorrZ-connValsZ);

% % Rank and rescale elements 
% maxVal = 1; minVal = -1;
% medPairCorrRank      = tiedrank(medPairCorr);
% medPairCorrRank      = minVal + ((medPairCorrRank-1)*(maxVal-minVal)/(length(medPairCorrRank)-1));

% % Calculate deviation from FC
% devPairCorr = abs(medPairCorr-connVals);
% devEnvelope = abs(medEnvelopeCorr-connVals);
% devInfraSlow = abs(medInfraSlowCorr-connVals);

% Calculate the ratio between FC and pairwise
% devPairCorr = abs(medPairCorr-connVals)./connVals;
% devEnvelope = abs(medEnvelopeCorr-connVals)./connVals;
% devInfraSlow = abs(medInfraSlowCorr-connVals)./connVals;

angleVals = abs(atand(connVals./medPairCorr) - 45);

figure;boxplot(angleVals,bandLabels);
ylabel('Angle (degrees)'); box off; 
[p,t,stats] = anova1(angleVals,bandLabels);
[c,m,h,gnames] = multcompare(stats);

%% Plot the distribution of deviation
figure; 
for iPlot = 1:3
    switch iPlot
        case 1
            plotVal = devPairCorr;
            plotLbl = 'Time series';
        case 2
            plotVal = devEnvelope;
            plotLbl = 'Power';
        case 3
            plotVal = devInfraSlow;
            plotLbl = 'Infra-slow';
    end
    subplot(1,3,iPlot);
    boxplot(plotVal,bandLabels); 
    ylabel('Bias');
    % violin(plotVal);
    title(plotLbl); box off; ylim([-0.05 1]);
end

[p,t,stats] = anova1(devPairCorr,bandLabels);
[c,m,h,gnames] = multcompare(stats);

%% Calculate the concordance correlation coefficient
for iP = 1:3
    switch iP
        case 1
            val = medPairCorr;
        case 2
            val = medEnvelopeCorr;
        case 3
            val = medInfraSlowCorr;
    end
    % 
    % r = mean((val-mean(val,1)).*(connVals-mean(connVals)),1,'omitnan') ./ sqrt(var(val,1)*var(connVals,1)); 
    % v = sqrt(var(val,1) ./ var(connVals,1));
    % u = (mean(val,1) - mean(connVals)) ./ sqrt(sqrt(var(val,1)) .* sqrt(var(connVals,1)));
    % biasFactor = 2 ./ (v + 1./v + u.^2);
    % ccc = r .* biasFactor;

        vif(:,iP) = diag(inv(corrcoef(vals)));

    y = zscore(connVals);
    X  = zscore(vals);
    clear  beta sigma eVal covMat

    % Get the linear model 
    mdlAll{iP} = fitlm(X,y,'VarNames',{'Theta','Alpha','Beta','Gamma','Spk','FC'});

    % Relative weights analysis to determine contributions of predictors
    [relImp(iP), r2(iP)] = rwa(X,y,bandLabels');

    % Dominance analysis
    [relativeImportance(:,iP),rsqDominance(iP)] = dominance(X,y);
    percentImportance(:,iP) = 100*relativeImportance(:,iP)./rsqDominance(iP);

end

% for iType = 1:3
%     switch iType
%         case 1
%             vals = medPairCorr;
%         case 2
%             vals = medEnvelopeCorr;
%         case 3
%             vals = medInfraSlowCorr;
%     end
% 
%     % Check for multicollinearity - uncomment if you want to visualize the
%     % correlations between frequencies 
%     % [~,~,h] = corrplot(double(vals),'VarNames',bandLabels,Type='Pearson',TestR="on");
%     % hAxes = findobj('Type','axes');
%     % set(hAxes(1:30),"XLim",[-0.6 1],"YLim",[-0.6 1]);
%     % cla(hAxes([1:5 6:6:30]));
% 
% 
%     % Checking for multicollinearity by calculating Variance Inflation
%     % factor. VIF = diagonal elements of inverse of correlation matrix.
%     vif(:,iType) = diag(inv(corrcoef(vals)));
% 
%     % vals = (vals-mean(vals))./std(vals);
%     y = zscore(connVals);
%     X  = zscore([vals distVals]);
%     clear  beta sigma eVal covMat
% 
%     % Get the linear model 
%     % mdlAll{iType} = fitlm(X,y,'VarNames',{'Theta','Alpha','Beta','Gamma','Spk','FC'});
% 
%     % Relative weights analysis to determine contributions of predictors
%     [relImp(iType), r2(iType)] = rwa(X,y,bandLabels');
% 
%     % Dominance analysis
%     [relativeImportance(:,iType),rsqDominance(iType)] = dominance(X,y);
%     percentImportance(:,iType) = 100*relativeImportance(:,iType)./rsqDominance(iType);
% end

% PCA for theta-beta frequencies
Z = zscore(medEnvelopeCorr); 
y = zscore(connVals);
[coeff,score,~,~,explained] = pca(Z(:,1:3));
pc1 = score(:,1); 
mdl = fitlm([pc1 Z(:,4:5)], zscore(connVals));
% [~,se,coeff1] = hac(mdl,'type','HC','weights','HC3');

% % Ridge regression
% lambdas = 0:0.1:20; B = ridge(y, Z, lambdas, 0); 
% figure; plot(lambdas, B);

resid_2 = Z(:,2) - [ones(size(Z,1),1) Z(:,3)]*regress(Z(:,2),[ones(size(Z,1),1) Z(:,3)]); 
resid_1 = Z(:,1) - [ones(size(Z,1),1) Z(:,3) Z(:,2)]*regress(Z(:,1),[ones(size(Z,1),1) Z(:,3) Z(:,2)]); 
mdl = fitlm([resid_1 resid_2 Z(:,3:5)], y);

%%

lambdas = 0:0.5:50;
cvp = cvpartition(length(y),'KFold',10);
mse = zeros(size(lambdas));
for i = 1:length(lambdas)
    errs = zeros(10,1);
    for j = 1:10
        tr = training(cvp,j); te = test(cvp,j);
        b = ridge(y(tr), [Z(tr,:)], lambdas(i), 0);
        pred = [ones(sum(te),1) Z(te,:)] * b;
        errs(j) = mean((y(te) - pred).^2);
    end
    mse(i) = mean(errs);
end
[~, idx] = min(mse); 
best_lambda = lambdas(idx);
b_final = ridge(y, Z, best_lambda, 0);


yhat = [ones(length(y),1) Z] * b_final; 
SS_res = sum((y - yhat).^2); SS_tot = sum((y - mean(y)).^2); 
R2_ridge = 1 - SS_res/SS_tot;

%% Pairwise correlations vs FC 
showScatterPlots(connVals,medPairCorr,medEnvelopeCorr,medInfraSlowCorr,...
    'Functional Connectivity','Pairwise correlations',[-0.6 1],[-0.1 1],...
    -0.5,0.9,0.8,bandLabels); 

%% Pairwise correlations vs distance
showScatterPlots(distVals,medPairCorr,medEnvelopeCorr,medInfraSlowCorr,'Distance',...
    'Pairwise correlations',[0 20],[-0.1 1],12,0.9,0.8,bandLabels);


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
    for iBand = 1:6
        [rhoVal(iType,iBand),pVal(iType,iBand)] = corr(connVals,var(:,iBand),"Type","Spearman");
    end
end

figure; 
imagesc(rhoVal); colorbar; clim([0 0.5]);
xticks(1:6); yticks(1:3);
xticklabels(bandLabels); yticklabels(timeLabels);

%% Plot 
for iType = 1:3
    switch iType
        case 1
            vals = medPairCorr;
        case 2
            vals = medEnvelopeCorr;
        case 3
            vals = medInfraSlowCorr;
    end

       % Checking for multicollinearity by calculating Variance Inflation
    % factor. VIF = diagonal elements of inverse of correlation matrix.
    vif(:,iType) = diag(inv(corrcoef(vals(:,1:5))));
    
    % vals = (vals-mean(vals))./std(vals);
    y = zscore(connVals);
    X  = zscore(vals);
    clear  beta sigma eVal covMat

    % Get the linear model 
    mdlAll{iType} = fitlm(X(:,1:5),y,'VarNames',{'Theta','Alpha','Beta','Gamma','Spk','FC'});
    
    % Relative weights analysis to determine contributions of predictors
    [relImp(iType), r2(iType)] = rwa(X,y,bandLabels(1:5)');

    % Dominance analysis
    [relativeImportance(:,iType),rsqDominance(iType)] = dominance(X,y);
    percentImportance(:,iType) = 100*relativeImportance(:,iType)./rsqDominance(iType);

    % Considering only 0-30 Hz, gamma, and spiking bands
    mdlThree{iType} = fitlm(X(:,4:6),y,'VarNames',{'Gamma','Spk','0-30 Hz','FC'});
    [relImpThree(iType), r2Three(iType)] = rwa(X(:,3:5),y,bandLabels(3:6)');

    % Dominance analysis
    [relativeImportanceThree(:,iType),rsqDominanceThree(iType)] = dominance(X(:,4:6),y);
    percentImportanceThree(:,iType) = 100*relativeImportanceThree(:,iType)./rsqDominanceThree(iType);

    vifThree(:,iType) = diag(inv(corrcoef(vals(:,4:6))));
end


%% Power across all channels within an electrode
meanSpecA = cat(1,allMonkeyVars(:).meanSpecValsAR); meanSpecA = meanSpecA./max(meanSpecA,[],3);
meanSpecB = cat(1,allMonkeyVars(:).meanSpecValsBR);meanSpecB = meanSpecB./max(meanSpecB,[],3);
meanSpecAll = (meanSpecA+meanSpecB)./2;

idx =1;
figure;
for iElec = 1:2

    switch iElec
        case 1
            vars = meanSpecA;
        case 2
            vars = meanSpecB;
    end
    medVal = squeeze(median(vars,1,'omitnan'));  % 5x21
    sem    = squeeze(mad(vars,1))./sqrt(size(vars,1));  % 5x21

    % Plot each frequency
    for iFreq = 1:5
        subplot(2,5,idx);
        medVal_freq = medVal(iFreq, :);  % 1x21
        sem_freq = sem(iFreq, :);        % 1x21

        plot(medVal_freq, 1:21); hold on;
        patch([medVal_freq-2*sem_freq, fliplr(medVal_freq+2*sem_freq)], ...
            [(1:21) fliplr(1:21)],'blue','FaceAlpha', 0.3, 'EdgeColor', 'none');
        idx = idx+1; title(bandLabels{iFreq});xlim([0.2 1]);
        box off; set(gca,'YDir','reverse');
        xlabel('Power (au)'); ylabel('Channels');
    end

end

%% Split power into superficial, middle,and deep compartments
figure;
for iBand = 1:5
    subplot(2,3,iBand); 
    
      yData = [median(squeeze(meanSpecAll(:,iBand,1:6)),2,'omitnan')...
        median(squeeze(meanSpecAll(:,iBand,7:12)),2,'omitnan') ...
        median(squeeze(meanSpecAll(:,iBand,13:end)),2,'omitnan')];
    boxplot(yData);ylabel('Power (a.u)')
    
    title(bandLabels{iBand}); ylim([0 1]); 
    xticklabels({'S','M','D'}); axis square; box off;
end

%% Compiling laminar data
% Finding recordings with single channels
singleChRow = [cellfun(@(x) isscalar(x),allMonkeyVars(1).intraCorrBR(:,1));...
    cellfun(@(x) isscalar(x),allMonkeyVars(2).intraCorrBR(:,1))];

% Time series
superAllPairCorr    = [allMonkeyVars(1).medPairCorrSuperR; allMonkeyVars(2).medPairCorrSuperR];
superAllPairCorr(singleChRow,:) = [];

midAllPairCorr  = [allMonkeyVars(1).medPairCorrMidR; allMonkeyVars(2).medPairCorrMidR];
midAllPairCorr(singleChRow,:) = [];

deepAllPairCorr = [allMonkeyVars(1).medPairCorrDeepR; allMonkeyVars(2).medPairCorrDeepR];
deepAllPairCorr(singleChRow) = []; 

superMidPair = [(allMonkeyVars(1).medPairASuperBMidR + allMonkeyVars(1).medPairAMidBSuperR)./2 ; ...
    (allMonkeyVars(2).medPairASuperBMidR + allMonkeyVars(2).medPairAMidBSuperR)./2];
superMidPair(singleChRow,:) = []; % Removing single channel recordings

superDeepPair = [(allMonkeyVars(1).medPairASuperBDeepR + allMonkeyVars(1).medPairADeepBSuperR)./2 ; ...
    (allMonkeyVars(2).medPairASuperBDeepR + allMonkeyVars(2).medPairADeepBSuperR)./2];
superDeepPair(singleChRow,:) = [];

midDeepPair = [(allMonkeyVars(1).medPairAMidBDeepR + allMonkeyVars(1).medPairADeepBMidR)./2 ; ...
    (allMonkeyVars(2).medPairAMidBDeepR + allMonkeyVars(2).medPairADeepBMidR)./2];
midDeepPair(singleChRow,:) = [];

% Power
superAllPow = [allMonkeyVars(1).envelopePairCorrSuperR; allMonkeyVars(2).envelopePairCorrSuperR];
superAllPow(singleChRow,:) = [];

midAllPow  = [allMonkeyVars(1).envelopePairCorrMidR; allMonkeyVars(2).envelopePairCorrMidR];
midAllPow(singleChRow,:) = [];

deepAllPow  = [allMonkeyVars(1).envelopePairCorrDeepR; allMonkeyVars(2).envelopePairCorrDeepR];
deepAllPow(singleChRow,:) = [];

superMidPairPow = [(allMonkeyVars(1).envelopeASuperBMidR + allMonkeyVars(1).envelopeAMidBSuperR)./2 ; ...
    (allMonkeyVars(2).envelopeASuperBMidR + allMonkeyVars(2).envelopeAMidBSuperR)./2];
superMidPairPow(singleChRow,:) = [];

superDeepPairPow = [(allMonkeyVars(1).envelopeADeepBSuperR + allMonkeyVars(1).envelopeASuperBDeepR)./2 ; ...
    (allMonkeyVars(2).envelopeADeepBSuperR + allMonkeyVars(2).envelopeASuperBDeepR)./2];
superDeepPairPow(singleChRow,:) = [];

midDeepPairPow = [(allMonkeyVars(1).envelopeAMidBDeepR + allMonkeyVars(1).envelopeADeepBMidR)./2 ; ...
    (allMonkeyVars(2).envelopeAMidBDeepR + allMonkeyVars(2).envelopeADeepBMidR)./2];
midDeepPairPow(singleChRow,:) = [];

% Infra-slow power 
superAllInfraPow = [allMonkeyVars(1).infraPairCorrSuperR; allMonkeyVars(2).infraPairCorrSuperR];
superAllInfraPow(singleChRow,:) = [];

midAllInfraPow  = [allMonkeyVars(1).infraPairCorrMidR; allMonkeyVars(2).infraPairCorrMidR];
midAllInfraPow(singleChRow,:) = [];

deepAllInfraPow  = [allMonkeyVars(1).infraPairCorrDeepR; allMonkeyVars(2).infraPairCorrDeepR];
deepAllInfraPow(singleChRow,:) = [];

superMidPairInfra = [(allMonkeyVars(1).infraASuperBMidR + allMonkeyVars(1). infraAMidBSuperR)./2 ; ...
    (allMonkeyVars(2).infraASuperBMidR + allMonkeyVars(2). infraAMidBSuperR)./2];
superMidPairInfra(singleChRow,:) = [];

superDeepPairInfra = [(allMonkeyVars(1).infraASuperBDeepR + allMonkeyVars(1).infraADeepBSuperR)./2 ; ...
    (allMonkeyVars(2).infraASuperBDeepR + allMonkeyVars(2).infraADeepBSuperR)./2];
superDeepPairInfra(singleChRow,:) = [];

midDeepPairInfra = [(allMonkeyVars(1).infraAMidBDeepR + allMonkeyVars(1).infraADeepBMidR)./2 ; ...
    (allMonkeyVars(2).infraAMidBDeepR + allMonkeyVars(2). infraADeepBMidR)./2];
midDeepPairInfra(singleChRow,:) = [];

% Removing single channel recordings from FC and distance
connValsNew = connVals; connValsNew(singleChRow)=[];
distValsNew = distVals; distValsNew(singleChRow) = [];

%% Plotting correspondence between pairwise correlations and FC 
% Looking at all types of correlations separately
laminarMat = [superAllPow superMidPairPow superDeepPairPow...
    superMidPairPow midAllPow midDeepPairPow ...
    superDeepPairPow midDeepPairPow deepAllPow ...
      ];

laminarMedMat = reshape(mean(laminarMat,1,'omitnan'),[6 3 3]);

figure; 
for iPlot = 1:6
    subplot(2,3,iPlot);
    imagesc(squeeze(laminarMedMat(iPlot,:,:))); colorbar; axis square;
    title(['Band: ' bandLabels{iPlot}]); 
    xticks(1:3); yticks(1:3); clim([0.1 0.4]);
    xticklabels(layerLabels); yticklabels(layerLabels);
end

% Organizing the data by compartment
laminarMargMat = [[superAllPow; superMidPairPow; superDeepPairPow]...
    [superMidPairPow; midAllPow; midDeepPairPow] ...
    [superDeepPairPow; midDeepPairPow; deepAllPow] ...
      ];

figure; 
for iPlot = 1:6
    subplot(2,3,iPlot);
    % boxplot((laminarMargMat(:,[iPlot iPlot+5 iPlot+10])),'Labels',{'S','M','D'}); 
    violin(laminarMargMat(:,[iPlot iPlot+6 iPlot+12]),...
        'xlabel',{'Superficial','Middle','Deep'},'bw',0.02,'edgecolor','none');
    box off;title(bandLabels{iPlot});
    xticklabels(layerLabels); %yticklabels(layerLabels);
    ylabel('Correlation'); legend off; ylim([-0.1 1]); yticks(0:0.2:1);
end


% Organizing and plotting the data for all frequencies and comparments in
% the same plot
laminarAllFreq = [reshape([superAllPow; superMidPairPow; superDeepPairPow],[],1) ...
        reshape([midAllPow; superMidPairPow; midDeepPairPow],[],1)...
        reshape([deepAllPow ; superDeepPairPow; midDeepPairPow],[],1)];

figure;violin(laminarAllFreq,...
    'xlabel',{'Superficial','Middle','Deep'},'bw',0.02,'edgecolor','none');

% Correspondence between FC and pairwise correlations 
laminarCorr = reshape(corr(connValsNew,laminarMat,'Type','Spearman'),[5 3 3]);
figure; 
for iPlot = 1:5
    subplot(2,3,iPlot);
    imagesc(squeeze(laminarCorr(iPlot,:,:))); colorbar; axis square;
    xticks(1:3); yticks(1:3); clim([0.1 0.5]);
    xticklabels(layerLabels); yticklabels(layerLabels);
end

%% Cortical area differences in correlations

pairClass(singleChRow,:) = [];
smLoc(singleChRow)       = [];
ssLoc(singleChRow)       = [];
mmLoc(singleChRow)       = [];

ssCorr = mean(reshape(corr(connValsNew(ssLoc),laminarMat(ssLoc,:),'Type','Spearman'),[5 3 3]),3,'omitnan');
mmCorr = mean(reshape(corr(connValsNew(mmLoc),laminarMat(mmLoc,:),'Type','Spearman'),[5 3 3]),3,'omitnan');
smCorr = mean(reshape(corr(connValsNew(smLoc),laminarMat(smLoc,:),'Type','Spearman'),[5 3 3]),3,'omitnan');

% Correspondence to FC
figure; 
for iPlot = 1:5
    plotVar = [ssCorr(iPlot,:); mmCorr(iPlot,:) ;smCorr(iPlot,:)];
    subplot(2,3,iPlot); imagesc(plotVar); colorbar; axis square;
    title(['Band: ' bandLabels{iPlot}]);
    xticks(1:3); xticklabels(layerLabels);
    yticks(1:3); yticklabels(areaLabels);
end

% Correlations 
ssMedMat = mean(reshape(mean(laminarMat(ssLoc,:),1,'omitnan'),[5 3 3]),3,'omitnan');
mmMedMat = mean(reshape(mean(laminarMat(mmLoc,:),1,'omitnan'),[5 3 3]),3,'omitnan');
smMedMat = mean(reshape(mean(laminarMat(smLoc,:),1,'omitnan'),[5 3 3]),3,'omitnan');

figure; 
for iPlot = 1:5
 plotVar = [ssMedMat(iPlot,:); mmMedMat(iPlot,:) ;smMedMat(iPlot,:)];
    subplot(2,3,iPlot); imagesc(plotVar); colorbar; axis square;
    title(['Band: ' bandLabels{iPlot}]);
    xticks(1:3); xticklabels(layerLabels);
    yticks(1:3); yticklabels(areaLabels);
end 

%% Group the data across animals and also for an individual animal 
for iM = 1:3
    if iM<3

        pairClass = allMonkeyVars(iM).pairClass;
        singleChRow = [cellfun(@(x) isscalar(x),allMonkeyVars(iM).intraCorrBR(:,1))];

        % Between probes
        % Time series
        superAll      = allMonkeyVars(iM).medPairCorrSuperR;
        midAll        = allMonkeyVars(iM).medPairCorrMidR;
        deepAll       = allMonkeyVars(iM).medPairCorrDeepR;
        superMidPair  = (allMonkeyVars(iM).medPairASuperBMidR + allMonkeyVars(iM).medPairAMidBSuperR)./2;
        superDeepPair = (allMonkeyVars(iM).medPairASuperBDeepR + allMonkeyVars(iM).medPairADeepBSuperR)./2 ;
        midDeepPair   = (allMonkeyVars(iM).medPairAMidBDeepR + allMonkeyVars(iM).medPairADeepBMidR)./2 ;

        % Power
        superPow         = allMonkeyVars(iM).envelopePairCorrSuperR;
        midPow           = allMonkeyVars(iM).envelopePairCorrMidR;
        deepPow          = allMonkeyVars(iM).envelopePairCorrDeepR;
        superMidPairPow  = (allMonkeyVars(iM).envelopeASuperBMidR + allMonkeyVars(iM).envelopeAMidBSuperR)./2 ;
        superDeepPairPow = (allMonkeyVars(iM).envelopeADeepBSuperR + allMonkeyVars(iM).envelopeASuperBDeepR)./2 ;
        midDeepPairPow   = (allMonkeyVars(iM).envelopeAMidBDeepR + allMonkeyVars(iM).envelopeADeepBMidR)./2 ;

        % Infraslow
        superInfra         = allMonkeyVars(iM).infraPairCorrSuperR;
        midInfra           = allMonkeyVars(iM).infraPairCorrMidR;
        deepInfra          = allMonkeyVars(iM).infraPairCorrDeepR;
        superMidPairInfra  = (allMonkeyVars(iM).infraASuperBMidR + allMonkeyVars(iM). infraAMidBSuperR)./2 ;
        superDeepPairInfra = (allMonkeyVars(iM).infraASuperBDeepR + allMonkeyVars(iM).infraADeepBSuperR)./2 ;
        midDeepPairInfra   = (allMonkeyVars(iM).infraAMidBDeepR + allMonkeyVars(iM).infraADeepBMidR)./2 ;

        connValsR = allMonkeyVars(iM).connValsR;
        distValsR = allMonkeyVars(iM).distValsR;

        % Within probes
        intraA = allMonkeyVars(iM).intraCorrAR;
        intraB = allMonkeyVars(iM).intraCorrBR;

        envelopeIntraA = allMonkeyVars(iM).envelopeIntraAAllR;
        envelopeIntraB = allMonkeyVars(iM).envelopeIntraBAllR;

        infraIntraA = allMonkeyVars(iM).infraIntraAAllR;
        infraIntraB = allMonkeyVars(iM).infraIntraBAllR;

        smLoc = sum(pairClass=='SM',2)==2;
        ssLoc = sum(pairClass=='SS',2)==2;
        mmLoc = sum(pairClass=='MM',2)==2;

        termination3b    = allMonkeyVars(iM).termination3b;
        terminationArea4 = allMonkeyVars(iM).terminationArea4;

    else
        pairClass = [allMonkeyVars(1).pairClass; allMonkeyVars(2).pairClass];
        singleChRow = [cellfun(@(x) isscalar(x),allMonkeyVars(1).intraCorrBR(:,1));...
            cellfun(@(x) isscalar(x),allMonkeyVars(2).intraCorrBR(:,1))];

        % Time series
        superAll      = cell2mat({allMonkeyVars.medPairCorrSuperR}');
        midAll        = cell2mat({allMonkeyVars.medPairCorrMidR}');
        deepAll       = cell2mat({allMonkeyVars.medPairCorrDeepR}');
        superMidPair  = (cell2mat({allMonkeyVars.medPairASuperBMidR}') + cell2mat({allMonkeyVars.medPairAMidBSuperR}'))./2;
        superDeepPair = (cell2mat({allMonkeyVars.medPairASuperBDeepR}') + cell2mat({allMonkeyVars.medPairADeepBSuperR}'))./2;
        midDeepPair   = (cell2mat({allMonkeyVars.medPairAMidBDeepR}') + cell2mat({allMonkeyVars.medPairADeepBMidR}'))./2;

        % Power
        superPow         = cell2mat({allMonkeyVars.envelopePairCorrSuperR}');
        midPow           = cell2mat({allMonkeyVars.envelopePairCorrMidR}');
        deepPow          = cell2mat({allMonkeyVars.envelopePairCorrDeepR}');
        superMidPairPow  = (cell2mat({allMonkeyVars.envelopeASuperBMidR}') + cell2mat({allMonkeyVars.envelopeAMidBSuperR}'))./2;
        superDeepPairPow = (cell2mat({allMonkeyVars.envelopeADeepBSuperR}') + cell2mat({allMonkeyVars.envelopeASuperBDeepR}'))./2;
        midDeepPairPow   = (cell2mat({allMonkeyVars.envelopeAMidBDeepR}') + cell2mat({allMonkeyVars.envelopeADeepBMidR}'))./2;

        % Infraslow
        superInfra         = cell2mat({allMonkeyVars.infraPairCorrSuperR}');
        midInfra           = cell2mat({allMonkeyVars.infraPairCorrMidR}');
        deepInfra          = cell2mat({allMonkeyVars.infraPairCorrDeepR}');
        superMidPairInfra  = (cell2mat({allMonkeyVars.infraASuperBMidR}') + cell2mat({allMonkeyVars. infraAMidBSuperR}'))./2;
        superDeepPairInfra = (cell2mat({allMonkeyVars.infraASuperBDeepR}') + cell2mat({allMonkeyVars.infraADeepBSuperR}'))./2;
        midDeepPairInfra   = (cell2mat({allMonkeyVars.infraAMidBDeepR}') + cell2mat({allMonkeyVars.infraADeepBMidR}'))./2;

        connValsR = [allMonkeyVars(1).connValsR; allMonkeyVars(2).connValsR];
        distValsR = [allMonkeyVars(1).distValsR; allMonkeyVars(2).distValsR];

        % Within probes
        intraA = [allMonkeyVars(1).intraCorrAR; allMonkeyVars(2).intraCorrAR];
        intraB = [allMonkeyVars(1).intraCorrBR; allMonkeyVars(2).intraCorrBR];

        envelopeIntraA = [allMonkeyVars(1).envelopeIntraAAllR; allMonkeyVars(2).envelopeIntraAAllR];
        envelopeIntraB = [allMonkeyVars(1).envelopeIntraBAllR; allMonkeyVars(2).envelopeIntraBAllR];

        infraIntraA = [allMonkeyVars(1).infraIntraAAllR;allMonkeyVars(2).infraIntraAAllR];
        infraIntraB = [allMonkeyVars(1).infraIntraBAllR;allMonkeyVars(2).infraIntraBAllR];

        smLoc = sum(pairClass=='SM',2)==2;
        ssLoc = sum(pairClass=='SS',2)==2;
        mmLoc = sum(pairClass=='MM',2)==2;

        termination3b    = cell2mat({allMonkeyVars.termination3b}');
        terminationArea4 = cell2mat({allMonkeyVars.terminationArea4}');
    end

    pairClass(singleChRow,:)      = [];
    smLoc(singleChRow)            = [];
    ssLoc(singleChRow)            = [];
    mmLoc(singleChRow)            = [];
    connValsR(singleChRow)        = [];
    distValsR(singleChRow)        = [];
    termination3b(singleChRow)    = [];
    terminationArea4(singleChRow) = [];

    intraA(singleChRow,:) = [];
    intraB(singleChRow,:) = [];
    envelopeIntraA(singleChRow,:) = [];
    envelopeIntraB(singleChRow,:) = [];
    infraIntraA(singleChRow,:) = [];
    infraIntraB(singleChRow,:) = [];

    % Average the within probe correlations
    avgIntraA = cell2mat(cellfun(@(x) mean(tril(x,-1),"all",'omitnan'),intraA,'UniformOutput',0));
    avgIntraB = cell2mat(cellfun(@(x) mean(tril(x,-1),"all",'omitnan'),intraB,'UniformOutput',0));

    avgEnvA = cell2mat(cellfun(@(x) mean(tril(x,-1),"all",'omitnan'),envelopeIntraA,'UniformOutput',0));
    avgEnvB = cell2mat(cellfun(@(x) mean(tril(x,-1),"all",'omitnan'),envelopeIntraB,'UniformOutput',0));

    avgInfraA = cell2mat(cellfun(@(x) mean(tril(x,-1),"all",'omitnan'),infraIntraA,'UniformOutput',0));
    avgInfraB = cell2mat(cellfun(@(x) mean(tril(x,-1),"all",'omitnan'),infraIntraB,'UniformOutput',0));

    % Plot the average correlations
    figure;
    subplot(131); boxplot([avgIntraA;avgIntraB],bandLabels); title('Within probe -time');box off; axis square; ylim([-0.05 0.5]);
    subplot(132); boxplot([avgEnvA; avgEnvB],bandLabels);    title('Within probe - power');box off; axis square; ylim([-0.05  0.5]);
    subplot(133); boxplot([avgInfraA; avgInfraB],bandLabels); title('Within probe - infraslow');box off; axis square; ylim([-0.05  0.5]);

    % Get the correlations corresponding to inter-laminar distance

    % intraDistA = cellfun(@(x) arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), ...
    %                                     1:size(x, 2)-1)', intraA, 'UniformOutput', false);
    chSize = 20;
    matSize = size(intraA);

    % Pairwise correlations in time
    intraDistA = cellfun(@(x) [arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), 1:(size(x, 2)-1)), ...
        NaN(1, chSize-(size(x, 2)-1))]',intraA,'UniformOutput',0);
    intraDistA = permute(reshape(permute(cat(3, intraDistA{:}), [3 1 2]), ...
        matSize(1), matSize(2), chSize), [1 3 2]);
    intraDistB = cellfun(@(x) [arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), 1:(size(x, 2)-1)), ...
        NaN(1, chSize-(size(x, 2)-1))]',intraB,'UniformOutput',0);
    intraDistB = permute(reshape(permute(cat(3, intraDistB{:}), [3 1 2]), ...
        matSize(1), matSize(2), chSize), [1 3 2]);

    % Power
    envDistA = cellfun(@(x) [arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), 1:(size(x, 2)-1)), ...
        NaN(1, chSize-(size(x, 2)-1))]',envelopeIntraA,'UniformOutput',0);
    envDistA = permute(reshape(permute(cat(3, envDistA{:}), [3 1 2]), ...
        matSize(1), matSize(2), chSize), [1 3 2]);
    envDistB = cellfun(@(x) [arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), 1:(size(x, 2)-1)), ...
        NaN(1, chSize-(size(x, 2)-1))]',envelopeIntraB,'UniformOutput',0);
    envDistB = permute(reshape(permute(cat(3, envDistB{:}), [3 1 2]), ...
        matSize(1), matSize(2), chSize), [1 3 2]);

    % Infra-slow power
    infraDistA = cellfun(@(x) [arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), 1:(size(x, 2)-1)), ...
        NaN(1, chSize-(size(x, 2)-1))]',infraIntraA,'UniformOutput',0);
    infraDistA = permute(reshape(permute(cat(3, infraDistA{:}), [3 1 2]), ...
        matSize(1), matSize(2), chSize), [1 3 2]);
    infraDistB = cellfun(@(x) [arrayfun(@(iX) mean(diag(x, -iX), 'omitnan'), 1:(size(x, 2)-1)), ...
        NaN(1, chSize-(size(x, 2)-1))]',infraIntraB,'UniformOutput',0);
    infraDistB = permute(reshape(permute(cat(3, infraDistB{:}), [3 1 2]), ...
        matSize(1), matSize(2), chSize), [1 3 2]);

    % Plot each type separately
    for iType = 1:3
        switch iType
            case 1
                vars = [intraDistA;intraDistB];
                typeLabel = 'Time series';
            case 2
                vars = [envDistA; envDistB];
                typeLabel = 'Power';
            case 3
                vars = [infraDistA; infraDistB];
                typeLabel = 'Infra-slow ';
        end

        figure;
        for iBand = 1:5
            subplot(2,3,iBand);
            % boxplot(squeeze(vars(:,:,iBand)));
            medVal = median(squeeze(vars(:,:,iBand)),1,'omitnan');
            sem    = mad(squeeze(vars(:,:,iBand)),1)./sqrt(size(vars,2));
            plot(1:20, medVal); hold on;
            patch([(1:20) fliplr(1:20)], [medVal-2*sem, fliplr(medVal+2*sem)], ...
                'blue', 'FaceAlpha', 0.3, 'EdgeColor', 'none');        title([bandLabels{iBand}]);
            xlabel('Distance (x100 um)'); box off; axis square;
            ylim([-0.2 1]);
            sgtitle(typeLabel);
        end
    end
  
    flagVal = true(size(connValsR)); %termination3b;% Change this to select specific pairs
    figure; idx = 1;
    for iRef = 1:3
        switch iRef
            case 1
                superVal = superAll;
                midVal   = midAll;
                deepVal  = deepAll;
                superMid = superMidPair;
                superDeep = superDeepPair;
                midDeep   = midDeepPair;
                typeLabel = 'Time series';

            case 2
                superVal = superPow ;
                midVal   = midPow ;
                deepVal  = deepPow ;
                superMid = superMidPairPow;
                superDeep = superDeepPairPow ;
                midDeep   = midDeepPairPow ;
                typeLabel = 'Power';

            case 3
                superVal = superInfra;
                midVal   = midInfra;
                deepVal  = deepInfra;
                superMid = superMidPairInfra;
                superDeep = superDeepPairInfra;
                midDeep   = midDeepPairInfra;
                typeLabel = 'Infraslow';
        end
        
        % Partial correlation between compartments and FC values
        [superAllCorr(iRef,:),pValSuperAll(iRef,:)] = partialcorr(superVal(flagVal,:),connValsR(flagVal),distValsR(flagVal),'Type','Spearman');
        [midAllCorr(iRef,:),pValMidAll(iRef,:)]     = partialcorr(midVal(flagVal,:),connValsR(flagVal),distValsR(flagVal),'Type','Spearman');
        [deepAllCorr(iRef,:),pValDeepAll(iRef,:)]   = partialcorr(deepVal(flagVal,:),connValsR(flagVal),distValsR(flagVal),'Type','Spearman');

        [superMidCorr(iRef,:),pValsSuperMid(iRef,:)]   = partialcorr(superMid(flagVal,:),connValsR(flagVal),distValsR(flagVal),'Type','Spearman');
        [superDeepCorr(iRef,:),pValsSuperDeep(iRef,:)] = partialcorr(superDeep(flagVal,:),connValsR(flagVal),distValsR(flagVal),'Type','Spearman');
        [midDeepCorr(iRef,:),pValsMidDeep(iRef,:)]     = partialcorr(midDeep(flagVal,:),connValsR(flagVal),distValsR(flagVal),'Type','Spearman');

        [superRefCorr(iRef,:),pSuperRef(iRef,:)] = partialcorr([superVal(flagVal,:); superMid(flagVal,:); superDeep(flagVal,:)],...
            repmat(connValsR(flagVal),[3 1]),repmat(distValsR(flagVal),[3 1]),'Type','Spearman');

        [midRefCorr(iRef,:),pMidRefCorr(iRef,:)] = partialcorr([midVal(flagVal,:); superMid(flagVal,:); midDeep(flagVal,:)],...
            repmat(connValsR(flagVal),[3 1]),repmat(distValsR(flagVal),[3 1]),'Type','Spearman');

        [deepRefCorr(iRef,:),pDeepRefCorr(iRef,:)] = partialcorr([deepVal(flagVal,:); midDeep(flagVal,:); superDeep(flagVal,:)],...
            repmat(connValsR(flagVal),[3 1]),repmat(distValsR(flagVal),[3 1]),'Type','Spearman');
        
        % Plot the correlations
        for iBand = 1:5
            subplot(3,5,idx); imagesc([median(superVal(flagVal,iBand)) median(superMid(flagVal,iBand)) median(superDeep(flagVal,iBand)); ...
                median(midVal(flagVal,iBand)) median(superMid(flagVal,iBand)) median(midDeep(flagVal,iBand)); ...
                median(deepVal(flagVal,iBand)) median(midDeep(flagVal,iBand)) median(superDeep(flagVal,iBand))]);%,'LineStyle','none');
            xticks(1:3); xticklabels({'Superficial','Middle','Deep'}); yticks(1:3); yticklabels({'Superficial','Middle','Deep'});
            axis square; colorbar; title([typeLabel '-' bandLabels{iBand}]); colormap parula;
            % maxVal = max([median(superVal(flagVal,iBand)) median(superMid(flagVal,iBand)) median(superDeep(flagVal,iBand)); ...
            %     median(midVal(flagVal,iBand)) median(superMid(flagVal,iBand)) median(midDeep(flagVal,iBand)); ...
            %     median(deepVal(flagVal,iBand)) median(midDeep(flagVal,iBand)) median(superDeep(flagVal,iBand))]);
            if iRef==2; clim([0 0.35]);
            elseif iRef ==3; clim([0 0.65]);
            else, clim([0 0.5]);
            end
            idx = idx+1;
        end
    end
end

%% Coherence for all frequencies
allCohRec = [allMonkeyVars.allCohRec]; % Coherence for all frequencies
allCerr   = {allMonkeyVars.allCerr}'; % Jackknife error intervals
allCerr   = cat(3,allCerr{:});
f         = allMonkeyVars(1).allCohVals.freqVals;
confCAvg  = allMonkeyVars(1).confCAvgBand;

% Plot coherence versus log frequency 
figure;
plot(f,median(allCohRec,2,'omitnan')); hold on;
  patch([f fliplr(f)], [median(squeeze(allCerr(1,:,:)),2,'omitnan');...
            flipud(median(squeeze(allCerr(2,:,:)),2,'omitnan'))],'blue','FaceAlpha',0.3,'EdgeColor','none')
  box off; xlabel(' Frequency (Hz)'); ylabel('Coherence');
  yline(confCAvg,'LineWidth',1); ylim([0 0.5]); set(gca, 'XScale', 'log');
 grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on');

% Organize the data into different frequency bands
clear cohValsFreq
connAll = connVals; 
distAll = distVals;
fs = 1e3;
thetaBand = [6 8];  alphaBand = [8 12]; betaBand  = [13 30]; gammaBand = [30 90];

for iBand = 1:4
    switch iBand
        case 1
            freqIdx = f>=thetaBand(1) & f<=thetaBand(2);
        case 2
            freqIdx = f>=alphaBand(1) & f<=alphaBand(2);
        case 3
            freqIdx = f>=betaBand(1) & f<=betaBand(2);
        case 4
            freqIdx = f>=gammaBand(1) & f<=gammaBand(2);
    end
    cohValsFreq(:,iBand) = median(allCohRec(freqIdx,:),1,'omitnan')';
end

% cohValsFreqR  = reshape(cohValsFreq,[maxRuns*maxDates 4]);

threshCond = sum(cohValsFreq<confCAvg,2)>=1;
cohValsFreq(threshCond,:) = []; connAll(threshCond) = []; distAll(threshCond) = [];

%% Plot coherence versus frequency
figure;
for iBand = 1:4
    subplot(2,2,iBand);
    plot(connAll,cohValsFreq(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
    ylim([0 0.8]); axis square;
    xVal = connAll; yVal = cohValsFreq(:,iBand);
    coeff = fit(xVal,double(yVal),'poly1','Robust','LAR');
    xFit  = linspace(min(xVal),max(xVal),1000);
    yFit  = coeff.p1*xFit + coeff.p2; mdl = fitlm(xVal,yVal,'RobustOpts','on');
    plot(xFit,yFit,'-k','LineWidth',1);
    text(0.8, 0.35,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(0.8, 0.3,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]); 
    title(bandLabels{iBand}); xlabel('Functional connectivity');
    ylabel('Coherence');
end

%% Plot coherence versus distance
figure;
for iBand = 1:4
    subplot(2,2,iBand);
    plot(distAll,cohValsFreq(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
    ylim([0 0.8]); axis square;
    xVal = distAll; yVal = cohValsFreq(:,iBand);

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


    text(8, 0.75,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(8, 0.7,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]);
    title(bandLabels{iBand}); xlabel('Distance (mm)');
    ylabel('Coherence');
end

y = zscore(connAll);
X  = zscore(cohValsFreq);
clear  beta sigma eVal covMat

% Get the linear model
mdlAll = fitlm(X,y,'VarNames',{'Theta','Alpha','Beta','Gamma','FC'});

% Dominance analysis: calculate the relative importance of each frequency
% band to the overall R squared value
[relativeImportance,rsqDominance] = dominance(X,y);
percentImportance = 100*relativeImportance./rsqDominance;

%% Phase-amplitude coupling
maxVal = 1; minVal = -1; 
% Load the modulation index values 
for iM = 1:2
   pacVals(iM) =  load(['D:\Data\' monkeys{iM} '_SqM\Left Hemisphere\allPACVars.mat']);
end

for iM = 1:3
    if iM<3 % Individual animal
        singleChRow = [cellfun(@(x) isscalar(x),allMonkeyVars(iM).intraCorrBR(:,1))];
        distValsR   = allMonkeyVars(iM).distValsR; distValsR(singleChRow) = [];
        connValsR   = allMonkeyVars(iM).connValsR; connValsR(singleChRow) = [];

        lowGammaThetaBtw  = pacVals(iM).thetaLGBtwT; % Low Gamma-theta between electrodes
        highGammaThetaBtw = pacVals(iM).thetaHGBtwT; % High Gamma-theta between electrodes

        lowGammaThetaWth  = pacVals(iM).thetaLGWthT; % Low Gamma-theta within electrodes
        highGammaThetaWth = pacVals(iM).thetaHGWthT; % High Gamma-theta within electrodes

        avgMIBtw = pacVals(iM).avgMIBetweenT; % Controls
        avgMIWth = pacVals(iM).avgMIWithinT; % Controls
        
        % Organized based on laminar compartments
        aAmpBPhaseLaminarT = pacVals.aAmpBPhaseLaminarT;
        aPhaseBAmpLaminarT = pacVals.aPhaseBAmpLaminarT;
        a2aLaminarT        = pacVals.a2aLaminarT;
        b2bLaminarT        = pacVals.b2bLaminarT;

        aAmpBPhaseThetaHighG = pacVals(iM).aAmpBPhaseThetaHighG;
        aAmpBPhaseThetaLowG  = pacVals(iM).aAmpBPhaseThetaLowG;

        aPhaseBAmpThetaHighG = pacVals(iM).aPhaseBAmpThetaHighG;
        aPhaseBAmpThetaLowG  = pacVals(iM).aPhaseBAmpThetaLowG;

        a2aLaminarThetaHighG = pacVals(iM).a2aLaminarThetaHighG;
        a2aLaminarThetaLowG  = pacVals(iM).a2aLaminarThetaLowG;

        b2bLaminarThetaHighG = pacVals(iM).b2bLaminarThetaHighG;
        b2bLaminarThetaLowG  = pacVals(iM).b2bLaminarThetaLowG;

        area3bARef = logical(allMonkeyVars(iM).termination3bARef);
        area3bBRef = logical(allMonkeyVars(iM).termination3bBRef);

        withinThresh = median([pacVals(iM).a2aThresh pacVals(iM).b2bThresh],'all');
        btwThresh    = median([pacVals(iM).a2bThresh pacVals(iM).b2aThresh],'all'); 
        
        termination3b = allMonkeyVars(iM).termination3b;

    else % Grouped data 
        singleChRow = [cellfun(@(x) isscalar(x),allMonkeyVars(1).intraCorrBR(:,1)); ...
           cellfun(@(x) isscalar(x),allMonkeyVars(2).intraCorrBR(:,1)) ];

        distValsR   = cell2mat({allMonkeyVars.distValsR}'); distValsR(singleChRow) = [];
        connValsR   = cell2mat({allMonkeyVars.connValsR}'); connValsR(singleChRow) = [];

        lowGammaThetaBtw  = cell2mat({pacVals.thetaLGBtwT}'); % Low Gamma-theta between electrodes
        highGammaThetaBtw = cell2mat({pacVals.thetaHGBtwT}'); % High Gamma-theta between electrodes

        lowGammaThetaWth  = cell2mat({pacVals.thetaLGWthT}'); % Low Gamma-theta within electrodes
        highGammaThetaWth = cell2mat({pacVals.thetaHGWthT}'); % High Gamma-theta within electrodes

        ctrlMIBtw  = cell2mat({pacVals.ctrlMIBtwT}'); % Controls
        ctrlMIWth  = cell2mat({pacVals.ctrlMIWithinT}'); % Controls

        avgMIBtw = cell2mat({pacVals.avgMIBetweenT}');
        avgMIWth = cell2mat({pacVals.avgMIWithinT}');

         % Organized based on laminar compartments
        aAmpBPhaseLaminarT = cell2mat({pacVals.aAmpBPhaseLaminarT}');
        aPhaseBAmpLaminarT = cell2mat({pacVals.aPhaseBAmpLaminarT}');
        a2aLaminarT        = cell2mat({pacVals.a2aLaminarT}');
        b2bLaminarT        = cell2mat({pacVals.b2bLaminarT}');

        aAmpBPhaseThetaHighG = cell2mat({pacVals.aAmpBPhaseThetaHighG}');
        aAmpBPhaseThetaLowG  = cell2mat({pacVals.aAmpBPhaseThetaLowG}');

        aPhaseBAmpThetaHighG = cell2mat({pacVals.aPhaseBAmpThetaHighG}');
        aPhaseBAmpThetaLowG  = cell2mat({pacVals.aPhaseBAmpThetaLowG}');

        a2aLaminarThetaHighG = cell2mat({pacVals.a2aLaminarThetaHighG}');
        a2aLaminarThetaLowG  = cell2mat({pacVals.a2aLaminarThetaLowG}');

        b2bLaminarThetaHighG = cell2mat({pacVals.b2bLaminarThetaHighG}');
        b2bLaminarThetaLowG  = cell2mat({pacVals.b2bLaminarThetaLowG}');

        area3bARef = logical(cell2mat({allMonkeyVars.termination3bARef}'));
        area3bBRef = logical(cell2mat({allMonkeyVars.termination3bBRef}'));

        termination3b    = cell2mat({allMonkeyVars.termination3b}');
        
        % Threshold
        withinThresh = median([pacVals.a2aThresh pacVals.b2bThresh],'all');
        btwThresh    = median([pacVals.a2bThresh pacVals.b2aThresh],'all');


    end
    termination3b(singleChRow) = [];

    % Plotting the modulation indices
    figure;
    subplot(121); boxplot([median(avgMIBtw,2,'omitnan') median(avgMIWth,2,'omitnan')],{'Between probes','Within probes'});
    box off; axis square; ylim([2.7e-4 4.3e-4]); ylabel('Average modulation index');
    subplot(122);violin([median(avgMIBtw,2,'omitnan') median(avgMIWth,2,'omitnan')],{'Between probes','Within probes'});
    box off; axis square; ylim([2.5e-4 4.5e-4]); ylabel('Average modulation index');


    figure;
    subplot(121); boxplot([lowGammaThetaBtw-btwThresh highGammaThetaBtw-btwThresh ...
        lowGammaThetaWth-withinThresh highGammaThetaWth-withinThresh],...
        {'Theta_LowGamma Between', 'Theta_highGamma Between','Theta_LowGamma Within',...
        'Theta_highGamma Within'});
    box off; axis square; ylim([2.4e-4 7e-4]);ylabel('Average modulation index');

    subplot(122); violin([lowGammaThetaBtw highGammaThetaBtw lowGammaThetaWth highGammaThetaWth],...
        {'Theta_LowGamma Between', 'Theta_highGamma Between','Theta_LowGamma Within',...
        'Theta_highGamma Within'});
    box off; axis square; ylim([2.4e-4 8e-4]);ylabel('Average modulation index');

    % Two way anova to test the effect of frequency pairs and MI
    % calculation
    rowNum = size(lowGammaThetaBtw,1);
    g1 = [repmat({'theta-lowG'},[rowNum 1]); repmat({'theta-highG'},[rowNum 1]); repmat({'theta-lowG'},[rowNum 1]);repmat({'theta-highG'},[rowNum 1])];
    g2 = [repmat({'Between'},[rowNum 1]); repmat({'Between'},[rowNum 1]); repmat({'Within'},[rowNum 1]);repmat({'Within'},[rowNum 1])];
    y  = reshape([[lowGammaThetaBtw; lowGammaThetaWth] [highGammaThetaBtw ;highGammaThetaWth]],[rowNum*4 1]);

   [p,~,stats,terms]= anovan(y,{g1,g2},'model','interaction','varnames',{'Frequency','Probes'});

   [results,~,~,gnames] = multcompare(stats,'Dimension',[1,2],'Alpha',0.001);
   tbl = array2table(results,"VariableNames", ...
       ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
   tbl.("Group A") = gnames(tbl.("Group A"));
   tbl.("Group B") = gnames(tbl.("Group B"));
 
   medMIBtw = median(avgMIBtw,2,'omitnan');
    
    idxMI = tiedrank(medMIBtw);
    idxMI = minVal + ((idxMI-1)*(maxVal-minVal)/(length(idxMI)-1));



    % lowGammaRank = (tiedrank(lowGammaThetaBtw)-0.5./(rowNum));
    % % % lowGammaRank = minVal + ((lowGammaRank-1)*(maxVal-minVal)/(length(lowGammaRank)-1));
    % highGammaRank = (tiedrank(highGammaThetaBtw)-0.5./(rowNum));
    % % % highGammaRank = minVal + ((highGammaRank-1)*(maxVal-minVal)/(length(highGammaRank)-1));
    % % 
    % ctrlRank      = (tiedrank(ctrlMIBtw)-0.5./(rowNum));
    % % % ctrlRank      = minVal + ((ctrlRank-1)*(maxVal-minVal)/(length(ctrlRank)-1));
    % % 
    % connValsRank = (tiedrank(connValsR)-0.5./(rowNum));
  epsilon = 0.1;
    lowGammaRank = (lowGammaThetaBtw-min(lowGammaThetaBtw))./(max(lowGammaThetaBtw)-min(lowGammaThetaBtw))* (1 - 2*epsilon) + epsilon;
    highGammaRank = (highGammaThetaBtw-min(highGammaThetaBtw))./(max(highGammaThetaBtw)-min(highGammaThetaBtw))* (1 - 2*epsilon) + epsilon;
    medPairCorrRank      = (ctrlMIBtw-min(ctrlMIBtw))./(max(ctrlMIBtw)-min(ctrlMIBtw))* (1 - 2*epsilon) + epsilon;
    connValsRank = (connValsR-min(connValsR))./(max(connValsR)-min(connValsR))* (1 - 2*epsilon) + epsilon;
    distValsRank = (distValsR-min(distValsR))./(max(distValsR)-min(distValsR))* (1-2*epsilon) + epsilon;

    figure;subplot(131);showLinearFit(connValsRank,lowGammaRank); axis square; title('Theta-low gamma');
    subplot(132); showLinearFit(connValsRank,highGammaRank); axis square; title('Theta-high gamma');
    subplot(133); showLinearFit(connValsRank,medPairCorrRank); axis square; title('Control');

    %%
    figure;subplot(131);showLinearFit(distValsRank,lowGammaRank); axis square; title('Theta-low gamma');
    subplot(132); showLinearFit(distValsRank,highGammaRank); axis square; title('Theta-high gamma');
    subplot(133); showLinearFit(distValsRank,medPairCorrRank); axis square; title('Control');

    
    %%
figure; 
 boxplot([[lowGammaThetaBtw(ssLoc) highGammaThetaBtw(ssLoc)] ...
   [lowGammaThetaBtw(mmLoc) highGammaThetaBtw(mmLoc); NaN(3,2) ]...
  [ lowGammaThetaBtw(smLoc) highGammaThetaBtw(smLoc); NaN(6,2)] ],...
   {'S1-S1:LowG','S1-S1: HighG','M1-M1:LowG','M1-M1: HighG','S1-M1:LowG','S1-M1: HighG'});

ylim([2.5e-4 6e-4]);box off;


%% 
figure; boxplot([ [lowGammaThetaBtw(ssLoc); highGammaThetaBtw(ssLoc)] ...
   [lowGammaThetaBtw(mmLoc);highGammaThetaBtw(mmLoc); NaN(6,1) ] ...
  [ lowGammaThetaBtw(smLoc);highGammaThetaBtw(smLoc); NaN(12,1)] ],...
  {'S1-S1','Motor-Motor','S1-Motor'});ylim([2.5e-4 6e-4]);box off;


    %%

    % figure;
    % miBtw = avgMIBtw;miBtw(isnan(miBtw))=0;
    % imagesc(imgaussian(reshape(median(miBtw,1),[nHigh nLow]),1));
    % % contourf(reshape(median(miBtw,1),[nHigh nLow]),'lines','none');
    %     set(gca,'YDir','normal'); clim([0 5e-4]);colorbar; axis square;shading interp
    %     % set(gca,'YDir','normal'); clim([0 6e-4]);colorbar; axis square;
    %     yticks(1:nHigh); xticks(1:nLow); yticklabels(gammaRange); xticklabels(lowFreqRange);
    %
    % figure;subplot(121); showLinearFit(connValsR,idxMI);axis square;
    % subplot(122);showLinearFit(connValsR,medMIBtw);axis square;
%%
vecSize = [length(connValsR) 3];

    % ampThetaLowG = [median([aAmpBPhaseThetaLowG(termination3b,1:3) ],2,'omitnan') ...
    %     median([aAmpBPhaseThetaLowG(termination3b,4:6)],2,'omitnan')...
    %     median([aAmpBPhaseThetaLowG(termination3b,7:9) ],2,'omitnan')];
    % % ampThetaLowG = tiedrank(reshape(ampThetaLowG,[numel(ampThetaLowG) 1]));
    % % ampThetaLowG = reshape(minVal + ((ampThetaLowG-1)*(maxVal-minVal)/(length(ampThetaLowG)-1)),vecSize);
    % 
    % phaseThetaLowG = [median([aAmpBPhaseThetaLowG(termination3b,[1 4 7]) ],2,'omitnan') ...
    %     median([aAmpBPhaseThetaLowG(termination3b,[2 5 8])] ,2,'omitnan')...
    %     median([aAmpBPhaseThetaLowG(termination3b,[3 6 9]) ],2,'omitnan')];
    %  % phaseThetaLowG = tiedrank(reshape(phaseThetaLowG,[numel(phaseThetaLowG) 1])); 
    %  % phaseThetaLowG = reshape(minVal + ((phaseThetaLowG-1)*(maxVal-minVal)/(length(phaseThetaLowG)-1)),vecSize);
    % 
    % ampThetaHighG = [median([aAmpBPhaseThetaHighG(termination3b,1:3) ],2,'omitnan') ...
    %     median([aAmpBPhaseThetaHighG(termination3b,4:6)],2,'omitnan')...
    %     median([aAmpBPhaseThetaHighG(termination3b,7:9) ],2,'omitnan')];
    %  % ampThetaHighG = tiedrank(reshape(ampThetaHighG,[numel(ampThetaHighG) 1])); 
    %  % ampThetaHighG = reshape(minVal + ((ampThetaHighG-1)*(maxVal-minVal)/(length(ampThetaHighG)-1)),vecSize);
    % 
    % phaseThetaHighG = [median([aAmpBPhaseThetaHighG(termination3b,[1 4 7])],2,'omitnan') ...
    %     median([aAmpBPhaseThetaHighG(termination3b,[2 5 8])],2,'omitnan')...
    %     median([aAmpBPhaseThetaHighG(termination3b,[3 6 9]) ...
    %     ],2,'omitnan')];


    ampThetaLowG = [median([aAmpBPhaseThetaLowG(area3bARef,1:3); aPhaseBAmpThetaLowG(area3bBRef,[1 4 7])],2,'omitnan') ...
        median([aAmpBPhaseThetaLowG(area3bARef,4:6) ;aPhaseBAmpThetaLowG(area3bBRef,[2 5 8])],2,'omitnan')...
        median([aAmpBPhaseThetaLowG(area3bARef,7:9); aPhaseBAmpThetaLowG(area3bBRef,[3 6 9])],2,'omitnan')];
    % ampThetaLowG = tiedrank(reshape(ampThetaLowG,[numel(ampThetaLowG) 1]));
    % ampThetaLowG = reshape(minVal + ((ampThetaLowG-1)*(maxVal-minVal)/(length(ampThetaLowG)-1)),vecSize);

    phaseThetaLowG = [median([aAmpBPhaseThetaLowG(area3bARef,[1 4 7]) ;aPhaseBAmpThetaLowG(area3bBRef,1:3)],2,'omitnan') ...
        median([aAmpBPhaseThetaLowG(area3bARef,[2 5 8]); aPhaseBAmpThetaLowG(area3bBRef,4:6)],2,'omitnan')...
        median([aAmpBPhaseThetaLowG(area3bARef,[3 6 9]); aPhaseBAmpThetaLowG(area3bBRef,7:9)],2,'omitnan')];
     % phaseThetaLowG = tiedrank(reshape(phaseThetaLowG,[numel(phaseThetaLowG) 1])); 
     % phaseThetaLowG = reshape(minVal + ((phaseThetaLowG-1)*(maxVal-minVal)/(length(phaseThetaLowG)-1)),vecSize);

    ampThetaHighG = [median([aAmpBPhaseThetaHighG(area3bARef,1:3) ;aPhaseBAmpThetaHighG(area3bBRef,[1 4 7])],2,'omitnan') ...
        median([aAmpBPhaseThetaHighG(area3bARef,4:6); aPhaseBAmpThetaHighG(area3bBRef,[2 5 8])],2,'omitnan')...
        median([aAmpBPhaseThetaHighG(area3bARef,7:9) ;aPhaseBAmpThetaHighG(area3bBRef,[3 6 9])],2,'omitnan')];
     % ampThetaHighG = tiedrank(reshape(ampThetaHighG,[numel(ampThetaHighG) 1])); 
     % ampThetaHighG = reshape(minVal + ((ampThetaHighG-1)*(maxVal-minVal)/(length(ampThetaHighG)-1)),vecSize);

    phaseThetaHighG = [median([aAmpBPhaseThetaHighG(area3bARef,[1 4 7]); aPhaseBAmpThetaHighG(area3bBRef,1:3)],2,'omitnan') ...
        median([aAmpBPhaseThetaHighG(area3bARef,[2 5 8]); aPhaseBAmpThetaHighG(area3bBRef,4:6)],2,'omitnan')...
        median([aAmpBPhaseThetaHighG(area3bARef,[3 6 9]) ;aPhaseBAmpThetaHighG(area3bBRef,7:9)],2,'omitnan')];
     % phaseThetaHighG = tiedrank(reshape(phaseThetaHighG,[numel(phaseThetaHighG) 1])); 
     % phaseThetaHighG = reshape(minVal + ((phaseThetaHighG-1)*(maxVal-minVal)/(length(phaseThetaHighG)-1)),vecSize);


    % ampThetaLowG = [median((aAmpBPhaseThetaLowG(:,1:3)+ aPhaseBAmpThetaLowG(:,[1 4 7]))/2,2,'omitnan') ...
    %     median((aAmpBPhaseThetaLowG(:,4:6)+aPhaseBAmpThetaLowG(:,[2 5 8]))/2,2,'omitnan')...
    %     median((aAmpBPhaseThetaLowG(:,7:9)+ aPhaseBAmpThetaLowG(:,[3 6 9]))/2,2,'omitnan')];
    % % ampThetaLowG = tiedrank(reshape(ampThetaLowG,[numel(ampThetaLowG) 1]));
    % % ampThetaLowG = reshape(minVal + ((ampThetaLowG-1)*(maxVal-minVal)/(length(ampThetaLowG)-1)),vecSize);
    % 
    % phaseThetaLowG = [median((aAmpBPhaseThetaLowG(:,[1 4 7])+aPhaseBAmpThetaLowG(:,1:3))/2,2,'omitnan') ...
    %     median((aAmpBPhaseThetaLowG(:,[2 5 8])+aPhaseBAmpThetaLowG(:,4:6))/2,2,'omitnan')...
    %     median((aAmpBPhaseThetaLowG(:,[3 6 9])+ aPhaseBAmpThetaLowG(:,7:9))/2,2,'omitnan')];
    %  % phaseThetaLowG = tiedrank(reshape(phaseThetaLowG,[numel(phaseThetaLowG) 1])); 
    %  % phaseThetaLowG = reshape(minVal + ((phaseThetaLowG-1)*(maxVal-minVal)/(length(phaseThetaLowG)-1)),vecSize);
    % 
    % ampThetaHighG = [median((aAmpBPhaseThetaHighG(:,1:3) +aPhaseBAmpThetaHighG(:,[1 4 7]))/2,2,'omitnan') ...
    %     median((aAmpBPhaseThetaHighG(:,4:6)+aPhaseBAmpThetaHighG(:,[2 5 8]))/2,2,'omitnan')...
    %     median((aAmpBPhaseThetaHighG(:,7:9)+aPhaseBAmpThetaHighG(:,[3 6 9]))/2,2,'omitnan')];
    %  % ampThetaHighG = tiedrank(reshape(ampThetaHighG,[numel(ampThetaHighG) 1])); 
    %  % ampThetaHighG = reshape(minVal + ((ampThetaHighG-1)*(maxVal-minVal)/(length(ampThetaHighG)-1)),vecSize);
    % 
    % phaseThetaHighG = [median((aAmpBPhaseThetaHighG(:,[1 4 7])+aPhaseBAmpThetaHighG(:,1:3))/2,2,'omitnan') ...
    %     median((aAmpBPhaseThetaHighG(:,[2 5 8])+aPhaseBAmpThetaHighG(:,4:6))/2,2,'omitnan')...
    %     median((aAmpBPhaseThetaHighG(:,[3 6 9])+aPhaseBAmpThetaHighG(:,7:9))/2,2,'omitnan')];
    %  phaseThetaHighG = tiedrank(reshape(phaseThetaHighG,[numel(phaseThetaHighG) 1])); 

%%
figure;
subplot(121); boxplot([ampThetaLowG phaseThetaLowG]);ylim([2.8e-4 7.2e-4]);%ylim([-1 1])%
title('Theta to low gamma'); box off; axis square; xticklabels(repmat({'Superficial','Middle','Deep'},1,2));

subplot(122); boxplot([ampThetaHighG phaseThetaHighG ]);ylim([3.5e-4 6.5e-4]);%ylim([-1 1]); %
title('Theta to high gamma'); box off; axis square; xticklabels(repmat({'Superficial','Middle','Deep'},1,2));

%%
figure; imagesc([corr(ampThetaLowG,connValsR(termination3b),'Type','Spearman') ...
    corr(phaseThetaLowG,connValsR(termination3b),'Type','Spearman') ...
    corr(ampThetaHighG,connValsR(termination3b),'Type','Spearman') ...
    corr(phaseThetaHighG,connValsR(termination3b),'Type','Spearman') ]);xticks(1:4);yticks(1:3);
yticklabels({'Superficial','Middle','Deep'}); colorbar
xticklabels({'Amplitude low G','Phase Low G','Amplitude High G','Phase High G'});

  %%
  reOrg = [1 4 7 2 5 8 3 6 9];
  distAmpRefLowG =[aAmpBPhaseThetaLowG(area3bARef,:); aPhaseBAmpThetaLowG(area3bBRef,reOrg)];
  distPhaseRefLowG = [aPhaseBAmpThetaLowG(area3bARef,:); aAmpBPhaseThetaLowG(area3bBRef,reOrg)];

  distAmpRefHighG = [aAmpBPhaseThetaHighG(area3bARef,:); aPhaseBAmpThetaHighG(area3bBRef,reOrg)];
  distPhaseRefHighG = [aPhaseBAmpThetaHighG(area3bARef,:); aAmpBPhaseThetaHighG(area3bBRef,reOrg)];

  ampRefLowG    = reshape(median([aAmpBPhaseThetaLowG(area3bBRef,:); aPhaseBAmpThetaLowG(area3bARef,reOrg)],1,'omitnan'),[3 3]);
  phaseRefLowG  = reshape(median([aPhaseBAmpThetaLowG(area3bBRef,:); aAmpBPhaseThetaLowG(area3bARef,reOrg)],1,'omitnan'),[3 3]);

  ampRefHighG    = reshape(median([aAmpBPhaseThetaHighG(area3bBRef,:); aPhaseBAmpThetaHighG(area3bARef,reOrg)],1,'omitnan'),[3 3]);
  phaseRefHighG  = reshape(median([aPhaseBAmpThetaHighG(area3bBRef,:); aAmpBPhaseThetaHighG(area3bARef,reOrg)],1,'omitnan'),[3 3]);

  % flagVal = ssLoc;
  %   distAmpRefLowG =[aAmpBPhaseThetaLowG(flagVal,:); aPhaseBAmpThetaLowG(flagVal,reOrg)];
  % distPhaseRefLowG = [aPhaseBAmpThetaLowG(flagVal,:); aAmpBPhaseThetaLowG(flagVal,reOrg)];
  % 
  % distAmpRefHighG = [aAmpBPhaseThetaHighG(flagVal,:); aPhaseBAmpThetaHighG(flagVal,reOrg)];
  % distPhaseRefHighG = [aPhaseBAmpThetaHighG(flagVal,:); aAmpBPhaseThetaHighG(flagVal,reOrg)];
  % 
  % 
  % ampRefLowG    = reshape(median([aAmpBPhaseThetaLowG(flagVal,:); aPhaseBAmpThetaLowG(flagVal,reOrg)],1,'omitnan'),[3 3]);
  % phaseRefLowG  = reshape(median([aPhaseBAmpThetaLowG(flagVal,:); aAmpBPhaseThetaLowG(flagVal,reOrg)],1,'omitnan'),[3 3]);
  % 
  % ampRefHighG    = reshape(median([aAmpBPhaseThetaHighG(flagVal,:); aPhaseBAmpThetaHighG(flagVal,reOrg)],1,'omitnan'),[3 3]);
  % phaseRefHighG  =  reshape(median([aPhaseBAmpThetaHighG(flagVal,:); aAmpBPhaseThetaHighG(flagVal,reOrg)],1,'omitnan'),[3 3]);

  figure;idx = 1;
    for iBand = 1:4
        switch iBand
            case 1
                plotVals = ampRefLowG;
                typeLabel = 'Amplitude - Low Gamma';
            case 2
                plotVals = ampRefHighG;
                typeLabel = 'Amplitude-High Gamma';
            case 3
                plotVals = phaseRefLowG;
                typeLabel = 'Phase-Low Gamma';
            case 4
                plotVals = phaseRefHighG;
                typeLabel = 'Phase-High Gamma';
        end

        subplot(2,2,idx);
        imagesc(plotVals);
        if iBand== 1|| iBand==3; clim([3.5e-4 4.5e-4]); else, clim([ 3.5e-4 5e-4]); end
        if iBand== 1|| iBand==3; clim([2.8e-4 5e-4]); else, clim([ 3.5e-4 5.5e-4]); end
        xticks(1:3); xticklabels({'Superficial','Middle','Deep'}); yticks(1:3); yticklabels({'Superficial','Middle','Deep'});
        colorbar; 
        axis square; colorbar; title(typeLabel); colormap parula;
        idx = idx+1;
    end
end

function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
% coeff = polyfit(xVal,yVal,1);
% xFit  = linspace(min(xVal),max(xVal),1000);
% yFit  = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
coeff = fit(xVal,double(yVal),'poly1','Robust','LAR');
xFit  = linspace(min(xVal),max(xVal),1000);
yFit  = coeff.p1*xFit + coeff.p2; mdl = fitlm(xVal,yVal,'RobustOpts','on');
plot(xFit,yFit,'-k','LineWidth',1);
if nargin<3
    textLocX  = max(xVal)-0.2*max(xVal);
    textLocY1 = max(yVal)-0.2*max(yVal);
    textLocY2 = max(yVal)-0.3*max(yVal);
end
text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end