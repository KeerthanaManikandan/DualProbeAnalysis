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

%% Combine all monkey data
for iM = 1:2
    allMonkeyVars{iM} =  load(['D:\Data\' monkeys{iM} '_SqM\Left Hemisphere\DualProbeVars.mat']);
end

connVals = [allMonkeyVars{1}.connValsR; allMonkeyVars{2}.connValsR];
distVals = [allMonkeyVars{1}.distValsR; allMonkeyVars{2}.distValsR];

% Pairwise correlations
medPairCorr      = [allMonkeyVars{1}.medPairCorrR; allMonkeyVars{2}.medPairCorrR];
medEnvelopeCorr  = [allMonkeyVars{1}.medCorrEnvelopeR; allMonkeyVars{2}.medCorrEnvelopeR];
medInfraSlowCorr = [allMonkeyVars{1}.medCorrInfraSlowR; allMonkeyVars{2}.medCorrInfraSlowR];

%%
bandLabels = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking'};
timeLabels = {'Time series','Power','Infraslow'};

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
mdl5 = fitlm(medEnvelopeCorr,connVals);

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

    % Check for multicollinearity - uncomment if you want to visualize the
    % correlations between frequencies 
    % [~,~,h] = corrplot(double(vals),'VarNames',bandLabels,Type='Pearson',TestR="on");
    % hAxes = findobj('Type','axes');
    % set(hAxes(1:30),"XLim",[-0.6 1],"YLim",[-0.6 1]);
    % cla(hAxes([1:5 6:6:30]));


    % Checking for multicollinearity by calculating Variance Inflation
    % factor. VIF = diagonal elements of inverse of correlation matrix.
    vif(:,iType) = diag(inv(corrcoef(vals)));
    
    % vals = (vals-mean(vals))./std(vals);
    y = zscore(connVals);
    X  = zscore(vals);
    clear  beta sigma eVal covMat

    % Get the linear model 
    mdlAll{iType} = fitlm(X,y,'VarNames',{'Theta','Alpha','Beta','Gamma','Spk','FC'});
    
    % Relative weights analysis to determine contributions of predictors
    [relImp(iType), r2(iType)] = rwa(X,y,bandLabels');

    % Dominance analysis
    [relativeImportance(:,iType),rsqDominance(iType)] = dominance(X,y);
    percentImportance(:,iType) = 100*relativeImportance(:,iType)./rsqDominance(iType);
end
%%
singleChRow = cellfun(@(x) isscalar(x),allMonkeyVars{2}.intraCorrBR(:,1));

superMidPair = [(allMonkeyVars{1}.medPairASuperBMidR + allMonkeyVars{1}.medPairAMidBSuperR)./2 ; ...
    (allMonkeyVars{2}.medPairASuperBMidR + allMonkeyVars{2}.medPairAMidBSuperR)./2];
superMidPair(find(singleChRow)+42,:) = [];

superDeepPair = [(allMonkeyVars{1}.medPairASuperBDeepR + allMonkeyVars{1}.medPairADeepBSuperR)./2 ; ...
    (allMonkeyVars{2}.medPairASuperBDeepR + allMonkeyVars{2}.medPairADeepBSuperR)./2];
superDeepPair(find(singleChRow)+42,:) = [];

midDeepPair = [(allMonkeyVars{1}.medPairAMidBDeepR + allMonkeyVars{1}.medPairADeepBMidR)./2 ; ...
    (allMonkeyVars{2}.medPairAMidBDeepR + allMonkeyVars{2}.medPairADeepBMidR)./2];
midDeepPair(find(singleChRow)+42,:) = [];

%
superMidPairPow = [(allMonkeyVars{1}.envelopeASuperBMidR + allMonkeyVars{1}.envelopeAMidBSuperR)./2 ; ...
    (allMonkeyVars{2}.envelopeASuperBMidR + allMonkeyVars{2}.envelopeAMidBSuperR)./2];
superMidPairPow(find(singleChRow)+42,:) = [];

superDeepPairPow = [(allMonkeyVars{1}.envelopeADeepBSuperR + allMonkeyVars{1}.envelopeASuperBDeepR)./2 ; ...
    (allMonkeyVars{2}.envelopeADeepBSuperR + allMonkeyVars{2}.envelopeASuperBDeepR)./2];
superDeepPairPow(find(singleChRow)+42,:) = [];

midDeepPairPow = [(allMonkeyVars{1}.envelopeAMidBDeepR + allMonkeyVars{1}.envelopeADeepBMidR)./2 ; ...
    (allMonkeyVars{2}.envelopeAMidBDeepR + allMonkeyVars{2}.envelopeADeepBMidR)./2];
midDeepPairPow(find(singleChRow)+42,:) = [];
%
superMidPairInfra = [(allMonkeyVars{1}.infraASuperBMidR + allMonkeyVars{1}. infraAMidBSuperR)./2 ; ...
    (allMonkeyVars{2}.infraASuperBMidR + allMonkeyVars{2}. infraAMidBSuperR)./2];
superMidPairInfra(find(singleChRow)+42,:) = [];

superDeepPairInfra = [(allMonkeyVars{1}.infraASuperBDeepR + allMonkeyVars{1}.infraADeepBSuperR)./2 ; ...
    (allMonkeyVars{2}.infraASuperBDeepR + allMonkeyVars{2}.infraADeepBSuperR)./2];
superDeepPairInfra(find(singleChRow)+42,:) = [];

midDeepPairInfra = [(allMonkeyVars{1}.infraAMidBDeepR + allMonkeyVars{1}.infraADeepBMidR)./2 ; ...
    (allMonkeyVars{2}.infraAMidBDeepR + allMonkeyVars{2}. infraADeepBMidR)./2];
midDeepPairInfra(find(singleChRow)+42,:) = [];


% figure;%('position',[932,103,709,838]);
connValsNew = connVals; connValsNew(find(singleChRow)+42)=[];
distValsNew = distVals; distValsNew(find(singleChRow)+42) = [];

corrSuperMid  = [corr(superMidPair,connValsNew) corr(superMidPairPow,connValsNew) corr(superMidPairInfra,connValsNew)];
corrSuperDeep = [corr(superDeepPair,connValsNew) corr(superDeepPairPow,connValsNew) corr(superDeepPairInfra,connValsNew)];
corrMidDeep   = [corr(midDeepPair,connValsNew) corr(midDeepPairPow,connValsNew) corr(midDeepPairInfra,connValsNew)];

% 
figure; % Organize by inter-laminar correlations
subplot(131); imagesc(corrSuperMid); title('Super-Middle');
xticks(1:3); xticklabels(timeLabels);yticks(1:5); yticklabels(bandLabels); axis square; clim([-0.2 0.5]); colormap jet; colorbar;

subplot(132); imagesc(corrSuperDeep);  title('Middle-Deep'); 
xticks(1:3); xticklabels(timeLabels);yticks(1:5); yticklabels(bandLabels); axis square; clim([-0.2 0.5]); colormap jet; colorbar;

subplot(133); imagesc(corrMidDeep);  title('Super-Deep');
xticks(1:3); xticklabels(timeLabels);yticks(1:5); yticklabels(bandLabels); axis square; clim([-0.2 0.5]); colormap jet; colorbar;

% Organize by timescale
figure; 
subplot(131); imagesc([corrSuperMid(:,1) corrMidDeep(:,1) corrSuperDeep(:,1)]'); title(timeLabels{1});
yticks(1:3); yticklabels({'S-M','M-D','S-D'});xticks(1:5); xticklabels(bandLabels); axis square; clim([-0.2 0.5]); colormap jet; colorbar;

subplot(132); imagesc([corrSuperMid(:,2) corrMidDeep(:,2) corrSuperDeep(:,2)]'); title(timeLabels{2});
yticks(1:3); yticklabels({'S-M','M-D','S-D'});xticks(1:5); xticklabels(bandLabels);axis square; clim([0 0.45]); colormap jet; colorbar;

subplot(133); imagesc([corrSuperMid(:,3) corrMidDeep(:,3) corrSuperDeep(:,3)]'); title(timeLabels{3});
yticks(1:3);yticklabels({'S-M','M-D','S-D'}); xticks(1:5); xticklabels(bandLabels);axis square; clim([-0.2 0.5]);colormap jet; colorbar;

% 
 figure; 
 subplot(131); boxplot([superMidPair(:,4) midDeepPair(:,4) superDeepPair(:,4)],{'S-M','M-D','S-D'});ylim([-0.2 1]);box off;
 subplot(132); boxplot([superMidPairPow(:,4) midDeepPairPow(:,4) superDeepPairPow(:,4)],{'S-M','M-D','S-D'});ylim([-0.2 1]);box off;
 subplot(133); boxplot([superMidPairInfra(:,4) midDeepPairInfra(:,4) superDeepPairInfra(:,4)],{'S-M','M-D','S-D'});ylim([-0.2 1]);box off;


% superMidPairCorrTheta = [ squeeze(allMonkeyVars{1}.superMedCorr(:,1,:))];
% superMidPairCorrAlpha = [squeeze(mVars.superMedCorr(:,1,:))];
% superMidPairCorrBeta  = [ squeeze(mVars.superMedCorr(:,2,:))];
% superMidPairCorrGamma = [ squeeze(mVars.superMedCorr(:,3,:))];
% superMidPairCorrSpk   = [ squeeze(mVars.superMedCorr(:,4,:))];

%     superMeanCorrEnvAlpha = [superMeanCorrEnvAlpha; squeeze(mVars.superMeanEnvelopeCorr(:,1,:))];
%     superMeanCorrEnvBeta  = [superMeanCorrEnvBeta; squeeze(mVars.superMeanEnvelopeCorr(:,2,:))];
%     superMeanCorrEnvGamma = [superMeanCorrEnvGamma; squeeze(mVars.superMeanEnvelopeCorr(:,3,:))];
%     superMeanCorrEnvWB    = [superMeanCorrEnvWB; squeeze(mVars.superMeanEnvelopeCorr(:,4,:))];
% 
%     superMedCorrEnvAlpha = [superMedCorrEnvAlpha; squeeze(mVars.superMedEnvelopeCorr(:,1,:))];
%     superMedCorrEnvBeta  = [superMedCorrEnvBeta; squeeze(mVars.superMedEnvelopeCorr(:,2,:))];
%     superMedCorrEnvGamma = [superMedCorrEnvGamma; squeeze(mVars.superMedEnvelopeCorr(:,3,:))];
%     superMedCorrEnvWB    = [superMedCorrEnvWB; squeeze(mVars.superMedEnvelopeCorr(:,4,:))];
% 
%     superMeanCorrInfraAlpha = [superMeanCorrInfraAlpha; squeeze(mVars.superMeanInfraSlowCorr(:,1,:))];
%     superMeanCorrInfraBeta  = [superMeanCorrInfraBeta; squeeze(mVars.superMeanInfraSlowCorr(:,2,:))];
%     superMeanCorrInfraGamma = [superMeanCorrInfraGamma; squeeze(mVars.superMeanInfraSlowCorr(:,3,:))];
%     superMeanCorrInfraWB    = [superMeanCorrInfraWB; squeeze(mVars.superMeanInfraSlowCorr(:,4,:))];
% 
%     superMedCorrInfraAlpha = [superMedCorrInfraAlpha; squeeze(mVars.superMedInfraSlowCorr(:,1,:))];
%     superMedCorrInfraBeta  = [superMedCorrInfraBeta; squeeze(mVars.superMedInfraSlowCorr(:,2,:))];
%     superMedCorrInfraGamma = [superMedCorrInfraGamma; squeeze(mVars.superMedInfraSlowCorr(:,3,:))];
%     superMedCorrInfraWB    = [superMedCorrInfraWB; squeeze(mVars.superMedInfraSlowCorr(:,4,:))];
% 
%     deepMeanPairCorrAlpha = [deepMeanPairCorrAlpha; squeeze(mVars.deepMeanCorr(:,1,:))];
%     deepMeanPairCorrBeta  = [deepMeanPairCorrBeta; squeeze(mVars.deepMeanCorr(:,2,:))];
%     deepMeanPairCorrGamma = [deepMeanPairCorrGamma; squeeze(mVars.deepMeanCorr(:,3,:))];
%     deepMeanPairCorrWB    = [deepMeanPairCorrWB; squeeze(mVars.deepMeanCorr(:,4,:))];
% 
%     deepMedPairCorrAlpha = [deepMedPairCorrAlpha; squeeze(mVars.deepMedCorr(:,1,:))];
%     deepMedPairCorrBeta  = [deepMedPairCorrBeta; squeeze(mVars.deepMedCorr(:,2,:))];
%     deepMedPairCorrGamma = [deepMedPairCorrGamma; squeeze(mVars.deepMedCorr(:,3,:))];
%     deepMedPairCorrWB    = [deepMedPairCorrWB; squeeze(mVars.deepMedCorr(:,4,:))];
% 
%     deepMeanCorrEnvAlpha = [deepMeanCorrEnvAlpha; squeeze(mVars.deepMeanEnvelopeCorr(:,1,:))];
%     deepMeanCorrEnvBeta  = [deepMeanCorrEnvBeta; squeeze(mVars.deepMeanEnvelopeCorr(:,2,:))];
%     deepMeanCorrEnvGamma = [deepMeanCorrEnvGamma; squeeze(mVars.deepMeanEnvelopeCorr(:,3,:))];
%     deepMeanCorrEnvWB    = [deepMeanCorrEnvWB; squeeze(mVars.deepMeanEnvelopeCorr(:,4,:))];
% 
%     deepMedCorrEnvAlpha = [deepMedCorrEnvAlpha; squeeze(mVars.deepMedEnvelopeCorr(:,1,:))];
%     deepMedCorrEnvBeta  = [deepMedCorrEnvBeta; squeeze(mVars.deepMedEnvelopeCorr(:,2,:))];
%     deepMedCorrEnvGamma = [deepMedCorrEnvGamma; squeeze(mVars.deepMedEnvelopeCorr(:,3,:))];
%     deepMedCorrEnvWB    = [deepMedCorrEnvWB; squeeze(mVars.deepMedEnvelopeCorr(:,4,:))];
% 
%     deepMeanCorrInfraAlpha = [deepMeanCorrInfraAlpha; squeeze(mVars.deepMeanInfraSlowCorr(:,1,:))];
%     deepMeanCorrInfraBeta  = [deepMeanCorrInfraBeta; squeeze(mVars.deepMeanInfraSlowCorr(:,2,:))];
%     deepMeanCorrInfraGamma = [deepMeanCorrInfraGamma; squeeze(mVars.deepMeanInfraSlowCorr(:,3,:))];
%     deepMeanCorrInfraWB    = [deepMeanCorrInfraWB; squeeze(mVars.deepMeanInfraSlowCorr(:,4,:))];
% 
%     deepMedCorrInfraAlpha = [deepMedCorrInfraAlpha; squeeze(mVars.deepMedInfraSlowCorr(:,1,:))];
%     deepMedCorrInfraBeta  = [deepMedCorrInfraBeta; squeeze(mVars.deepMedInfraSlowCorr(:,2,:))];
%     deepMedCorrInfraGamma = [deepMedCorrInfraGamma; squeeze(mVars.deepMedInfraSlowCorr(:,3,:))];
%     deepMedCorrInfraWB    = [deepMedCorrInfraWB; squeeze(mVars.deepMedInfraSlowCorr(:,4,:))];
% 


% for iM = 1:2
%     switch iM
%         case 1
%             monkeyName = 'CharlieSheen';
%         case 2
%             monkeyName = 'Whiskey';
%     end
% 
%     mVars     = load(['D:\Data\' monkeyName '_SqM\Left Hemisphere\DualProbeVars.mat']);
%     dist      = [dist; mVars.distValsR];
%     fcVals    = [fcVals; mVars.connValsR];
%     % pairType  = [pairType; mVars.pairClass];
%     % patchInfo = [patchInfo mVars.patchInfo];
% 
%     if iM == 1; lenBadDat = 0; else ; lenBadDat = size(meanPairCorrAlpha,1); end
%     % badDatIdx = [badDatIdx; mVars.nanIdx + lenBadDat];
%     % Dividing pairwise correlations based on superficial and deeper layers
%     % of cortex
%     for iList = 1:size(mVars.marginalA_B_temp,2)
%         if ~isempty(mVars.marginalA_B_temp{iList})
%             % Normalize the vector
%             numRows = min([size(mVars.marginalA_B_temp{iList},1) size(mVars.marginalB_A_temp{iList},1)]);
%             if numRows>=20
%                 numRows = 20;
%             else
%                 numRows = min([size(mVars.marginalA_B_temp{iList},1) size(mVars.marginalB_A_temp{iList},1)]);
%             end
%             alphaVarsAB(:,:,rowIdx) = squeeze([mVars.marginalA_B_temp{iList}(1:numRows,1,:)./max(mVars.marginalA_B_temp{iList}(1:numRows,1,:)) ;zeros(20-numRows,1,2)]);
%             alphaVarsBA(:,:,rowIdx) = squeeze([mVars.marginalB_A_temp{iList}(1:numRows,1,:)./max(mVars.marginalB_A_temp{iList}(1:numRows,1,:)) ;zeros(20-numRows,1,2)]);
% 
%             gammaVarsAB(:,:,rowIdx) = squeeze([mVars.marginalA_B_temp{iList}(1:numRows,3,:)./max(mVars.marginalA_B_temp{iList}(1:numRows,3,:)) ;zeros(20-numRows,1,2)]);
%             gammaVarsBA(:,:,rowIdx) = squeeze([mVars.marginalB_A_temp{iList}(1:numRows,3,:)./max(mVars.marginalB_A_temp{iList}(1:numRows,3,:)) ;zeros(20-numRows,1,2)]);
% 
%             one_oneMapAlpha(:,:,rowIdx) = squeeze([mVars.one_oneCorr_temp{iList}(1:numRows,1,:)./max(mVars.one_oneCorr_temp{iList}(1:numRows,1,:)) ;zeros(20-numRows,1,2)]);
%             one_oneMapGamma(:,:,rowIdx) = squeeze([mVars.one_oneCorr_temp{iList}(1:numRows,3,:)./max(mVars.one_oneCorr_temp{iList}(1:numRows,3,:)) ;zeros(20-numRows,1,2)]);
%             rowIdx = rowIdx+1;
%         else
%             removeIdx = [removeIdx rowIdx];
%             rowIdx = rowIdx+1;
% 
%         end
%     end
% end
% badDatIdx = [badDatIdx; 77;78;79]; % Remove anesthesia controls from data

%% 
superAllCorr = [allMonkeyVars{1}.envelopePairCorrSuperR; allMonkeyVars{2}.envelopePairCorrSuperR];
midAllCorr   = [allMonkeyVars{1}.envelopePairCorrMidR; allMonkeyVars{2}.envelopePairCorrMidR];
deepAllCorr  = [allMonkeyVars{1}.envelopePairCorrDeepR; allMonkeyVars{2}.envelopePairCorrDeepR];

superAllCorr(find(singleChRow)+42,:) = [];
midAllCorr(find(singleChRow)+42,:) = [];
deepAllCorr(find(singleChRow)+42,:) = [];


laminarCorr = [corr(superAllCorr,connValsNew) corr(midAllCorr,connValsNew) corr(deepAllCorr,connValsNew)];
figure; imagesc(laminarCorr'); colorbar
colormap jet;  axis square; clim([0 0.5]);


for iType = 1:3
    clear var
    switch iType
        case 1
            var = superAllCorr;
        case 2
            var = midAllCorr;
        case 3
            var = deepAllCorr; 
    end
    for iBand = 1:5
        [rhoValLaminar(iType,iBand),pVal(iType,iBand)] = corr(connValsNew,var(:,iBand),'Type','Spearman');
    end
end

figure; bar([laminarCorr(4,:) corrSuperMid(4,2) corrMidDeep(4,2) corrSuperDeep(4,2)]); 
box off; 

figure; boxplot([superAllCorr(:,4) midAllCorr(:,4) deepAllCorr(:,4) superMidPairPow(:,4) midDeepPairPow(:,4) superDeepPairPow(:,4)]);

%% 
clear corrSuperMid corrMidDeep corrSuperDeep laminarCorr
pairClass = [allMonkeyVars{1}.pairClass; allMonkeyVars{2}.pairClass];
pairClass(find(singleChRow)+42,:) =[];

smLoc = sum(pairClass=='SM',2)==2;
ssLoc = sum(pairClass=='SS',2)==2;
mmLoc = sum(pairClass=='MM',2)==2;

for iPlot = 1:3
    switch iPlot
        case 1
            val = ssLoc;
            figTitle ='Somatosensory-Somatosensory';
        case 2
            val = mmLoc;
            figTitle ='Motor-Motor';            
        case 3
            val = smLoc;
            figTitle ='Somatosensory-Motor';
    end
    allSuper = superAllCorr(val,:); 
    allMid   = midAllCorr(val,:);
    allDeep  = deepAllCorr(val,:);
    connNew  = connValsNew(val,:);
    distNew  = distValsNew(val,:);
    superMid = superMidPairPow(val,:);
    midDeep  = midDeepPairPow(val,:);
    superDeep = superDeepPairPow(val,:);
    
    [rAllSuper(iPlot,:),pAllSuper(iPlot,:)] = corr(allSuper,connNew,'Type','Spearman'); 
    [rAllMid(iPlot,:),pAllMid(iPlot,:)] = corr(allMid,connNew,'Type','Spearman');
    [rAllDeep(iPlot,:),pAllDeep(iPlot,:)] = corr(allDeep,connNew,'Type','Spearman');
    laminarCorr   = [corr(allSuper,connNew,'Type','Spearman') corr(allMid,connNew,'Type','Spearman') corr(allDeep,connNew,'Type','Spearman')]';
    [corrSuperMid(iPlot,:),pSuperMid(iPlot,:)]  = corr(superMid,connNew,'Type','Spearman') ;
    [corrMidDeep(iPlot,:),pMidDeep(iPlot,:)]   = corr(midDeep,connNew,'Type','Spearman');
    [corrSuperDeep(iPlot,:),pSuperDeep(iPlot,:)]  = corr(superDeep,connNew,'Type','Spearman');

    figure; 
    subplot(121);
    bar([laminarCorr(:,4)' corrSuperMid(iPlot,4) corrMidDeep(iPlot,4) corrSuperDeep(iPlot,4)]);
    box off; ylim([0 0.6])
    subplot(122); 
    boxplot([allSuper(:,4) allMid(:,4) allDeep(:,4) superMid(:,4) midDeep(:,4) superDeep(:,4)]);
    xticks(1:6); xticklabels({'Super-Super','Middle-Middle','Deep-Deep','Super-Mid','Mid-Deep','Super-Deep'});
   % swarmchart((1:6).*ones(size(allSuper,1),1),[allSuper(:,4) allMid(:,4) allDeep(:,4) superMid(:,4) midDeep(:,4) superDeep(:,4)],80,'Marker','.','YJitterWidth',0.2);
    box off; ylim([0 1]);
    sgtitle(figTitle);
    [p(iPlot,:),~,stats{iPlot}] = anova1([allSuper(:,4) allMid(:,4) allDeep(:,4) superMid(:,4) midDeep(:,4) superDeep(:,4)],[],'off');
    [c{iPlot}] = multcompare(stats{iPlot},'CriticalValueType','bonferroni','Display','off');
end

figure;imagesc([rAllSuper(1,4) rAllMid(1,4) rAllDeep(1,4) corrSuperMid(1,4) corrMidDeep(1,4) corrSuperDeep(1,4);...
    rAllSuper(2,4) rAllMid(2,4) rAllDeep(2,4) corrSuperMid(2,4) corrMidDeep(2,4) corrSuperDeep(2,4);...
    rAllSuper(3,4) rAllMid(3,4) rAllDeep(3,4) corrSuperMid(3,4) corrMidDeep(3,4) corrSuperDeep(3,4)]);
colormap turbo; colorbar; clim([0 0.5]); 
yticks(1:3); yticklabels({'S-S','M-M','S-M'}); 
xticks(1:6); xticklabels({'Super-Super','Middle-Middle','Deep-Deep','Super-Mid','Mid-Deep','Super-Deep'});