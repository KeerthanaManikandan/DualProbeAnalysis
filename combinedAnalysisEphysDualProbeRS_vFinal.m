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

monkeys     = {'CharlieSheen'; 'Whiskey'};
bandLabels  = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking'};
timeLabels  = {'Time series','Power','Infraslow'};
layerLabels = {'S','M','D'};
areaLabels  = {'S-S','M-M','S-M'};
hemisphere  = 'Left';

for iM = 1:2
    allMonkeyVars(iM) =  load(['D:\Data\' monkeys{iM} '_SqM\Left Hemisphere\DualProbeVars.mat']);
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

termination3b    = cell2mat({allMonkeyVars.termination3b}');
terminationArea4 = cell2mat({allMonkeyVars.terminationArea4}');


%% Pairwise correlations vs FC 
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

laminarMedMat = reshape(mean(laminarMat,1,'omitnan'),[5 3 3]);

figure; 
for iPlot = 1:5
    subplot(2,3,iPlot);
    imagesc(squeeze(laminarMedMat(iPlot,:,:))); colorbar; axis square;
    title(['Band: ' bandLabels{iPlot}]); 
    xticks(1:3); yticks(1:3); 
    xticklabels(layerLabels); yticklabels(layerLabels);
end

% Organizing the data by compartment
laminarMargMat = [[superAllPow; superMidPairPow; superDeepPairPow]...
    [superMidPairPow; midAllPow; midDeepPairPow] ...
    [superDeepPairPow; midDeepPairPow; deepAllPow] ...
      ];

figure; 
for iPlot = 1:5
    subplot(2,3,iPlot);
    % boxplot((laminarMargMat(:,[iPlot iPlot+5 iPlot+10])),'Labels',{'S','M','D'}); 
    violin(laminarMargMat(:,[iPlot iPlot+5 iPlot+10]),...
        'xlabel',{'Superficial','Middle','Deep'},'bw',0.02,'edgecolor','none');
    box off;title(bandLabels{iPlot});
    xticklabels(layerLabels); %yticklabels(layerLabels);
    ylabel('Correlation');
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
    xticks(1:3); yticks(1:3); 
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


%% Coherence

