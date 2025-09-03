%% combinedAnalysisEphysDualProbeRS
% This function performs analysis on LFP recorded simultaneously from two
% linear electrode arrays for TWO monkeys
% March 2, 2024 - Keerthana Manikandan
% Refer to pow_Infraslow.m for previous version

clear;clc;
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
saveFigureFlag = 0;
badChProbeB    = [14 22];
fs             = 1e3; % Sampling frequency
gammaBand      = [30 90]; [zG,pG,kG] = butter(3,gammaBand./(fs/2),'bandpass'); [sosG,gG] = zp2sos(zG,pG,kG);
alphaBand      = [8 12];  [zA,pA,kA] = butter(3,alphaBand./(fs/2),'bandpass'); [sosA,gA] = zp2sos(zA,pA,kA);
betaBand       = [13 30]; [zB,pB,kB] = butter(3,betaBand./(fs/2),'bandpass');  [sosB,gB] = zp2sos(zB,pB,kB);
% probeBList     = [1:13 15:21 23:32];
chOutCortex    = 1:3;
chDeep         = 30:32;

rowIdx = 1;
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

[z,p,k] = butter(2,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma';'Spiking'};

%% Load both monkey data
clear alphaVarsAB alphaVarsBA gammaVarsAB gammaVarsBA one_oneMapAlpha one_oneMapGamma
rowIdx = 1; dist = [];fcVals = []; patchInfo = []; pairType = []; badDatIdx = []; removeIdx = [];

meanPairCorrAlpha = []; medPairCorrAlpha =[]; meanCorrEnvAlpha =[]; medCorrEnvAlpha = []; meanCorrInfraAlpha = []; medCorrInfraAlpha = [];
superMeanPairCorrAlpha = []; superMedPairCorrAlpha = []; superMeanCorrEnvAlpha = []; superMedCorrEnvAlpha =[]; superMeanCorrInfraAlpha = []; superMedCorrInfraAlpha =[];
deepMeanPairCorrAlpha = [];  deepMedPairCorrAlpha = []; deepMeanCorrEnvAlpha = []; deepMedCorrEnvAlpha =[]; deepMeanCorrInfraAlpha = []; deepMedCorrInfraAlpha =[];

meanPairCorrBeta = []; medPairCorrBeta =[]; meanCorrEnvBeta =[]; medCorrEnvBeta = []; meanCorrInfraBeta = []; medCorrInfraBeta = [];
superMeanPairCorrBeta = []; superMedPairCorrBeta = []; superMeanCorrEnvBeta = []; superMedCorrEnvBeta =[]; superMeanCorrInfraBeta = []; superMedCorrInfraBeta =[];
deepMeanPairCorrBeta = [];  deepMedPairCorrBeta = []; deepMeanCorrEnvBeta = []; deepMedCorrEnvBeta =[]; deepMeanCorrInfraBeta = []; deepMedCorrInfraBeta =[];

meanPairCorrGamma = []; medPairCorrGamma =[]; meanCorrEnvGamma =[]; medCorrEnvGamma = []; meanCorrInfraGamma = []; medCorrInfraGamma = [];
superMeanPairCorrGamma = []; superMedPairCorrGamma = []; superMeanCorrEnvGamma = []; superMedCorrEnvGamma =[]; superMeanCorrInfraGamma = []; superMedCorrInfraGamma =[];
deepMeanPairCorrGamma = [];  deepMedPairCorrGamma = []; deepMeanCorrEnvGamma = []; deepMedCorrEnvGamma =[]; deepMeanCorrInfraGamma = []; deepMedCorrInfraGamma =[];

meanPairCorrWB = []; medPairCorrWB =[]; meanCorrEnvWB =[]; medCorrEnvWB = []; meanCorrInfraWB = []; medCorrInfraWB = [];
superMeanPairCorrWB = []; superMedPairCorrWB = []; superMeanCorrEnvWB = []; superMedCorrEnvWB =[]; superMeanCorrInfraWB = []; superMedCorrInfraWB =[];
deepMeanPairCorrWB = [];  deepMedPairCorrWB = []; deepMeanCorrEnvWB = []; deepMedCorrEnvWB =[]; deepMeanCorrInfraWB = []; deepMedCorrInfraWB =[];

for iM = 1:2
    switch iM
        case 1
            monkeyName = 'CharlieSheen';
        case 2
            monkeyName = 'Whiskey';
    end

    mVars     = load(['D:\Data\' monkeyName '_SqM\Left Hemisphere\DualProbeVars.mat']);
    dist      = [dist; mVars.distSitesAll];
    fcVals    = [fcVals; mVars.connValsAll];
    pairType  = [pairType; mVars.pairClass];
    patchInfo = [patchInfo mVars.patchInfo];

    if iM == 1; lenBadDat = 0; else ; lenBadDat = size(meanPairCorrAlpha,1); end
    badDatIdx = [badDatIdx; mVars.nanIdx + lenBadDat];

    % Pairwise correlations
    meanPairCorrAlpha = [meanPairCorrAlpha; squeeze(mVars.meanPairCorr(:,2,:))];
    meanPairCorrBeta  = [meanPairCorrBeta; squeeze(mVars.meanPairCorr(:,3,:))];
    meanPairCorrGamma = [meanPairCorrGamma; squeeze(mVars.meanPairCorr(:,4,:))];
    meanPairCorrWB    = [meanPairCorrWB; squeeze(mVars.meanPairCorr(:,4,:))];

    medPairCorrAlpha = [medPairCorrAlpha; squeeze(mVars.medPairCorr(:,2,:))];
    medPairCorrBeta  = [medPairCorrBeta; squeeze(mVars.medPairCorr(:,3,:))];
    medPairCorrGamma = [medPairCorrGamma; squeeze(mVars.medPairCorr(:,4,:))];
    medPairCorrWB    = [medPairCorrWB; squeeze(mVars.medPairCorr(:,5,:))];

    meanCorrEnvAlpha = [meanCorrEnvAlpha; squeeze(mVars.meanCorrEnvelope(:,1,:))];
    meanCorrEnvBeta  = [meanCorrEnvBeta; squeeze(mVars.meanCorrEnvelope(:,2,:))];
    meanCorrEnvGamma = [meanCorrEnvGamma; squeeze(mVars.meanCorrEnvelope(:,3,:))];
    meanCorrEnvWB    = [meanCorrEnvWB; squeeze(mVars.meanCorrEnvelope(:,4,:))];

    medCorrEnvAlpha = [medCorrEnvAlpha; squeeze(mVars.medCorrEnvelope(:,1,:))];
    medCorrEnvBeta  = [medCorrEnvBeta; squeeze(mVars.medCorrEnvelope(:,2,:))];
    medCorrEnvGamma = [medCorrEnvGamma; squeeze(mVars.medCorrEnvelope(:,3,:))];
    medCorrEnvWB    = [medCorrEnvWB; squeeze(mVars.medCorrEnvelope(:,4,:))];

    meanCorrInfraAlpha = [meanCorrInfraAlpha; squeeze(mVars.meanCorrInfraSlow(:,1,:))];
    meanCorrInfraBeta  = [meanCorrInfraBeta; squeeze(mVars.meanCorrInfraSlow(:,2,:))];
    meanCorrInfraGamma = [meanCorrInfraGamma; squeeze(mVars.meanCorrInfraSlow(:,3,:))];
    meanCorrInfraWB    = [meanCorrInfraWB; squeeze(mVars.meanCorrInfraSlow(:,4,:))];

    medCorrInfraAlpha = [medCorrInfraAlpha; squeeze(mVars.medCorrInfraSlow(:,1,:))];
    medCorrInfraBeta  = [medCorrInfraBeta; squeeze(mVars.medCorrInfraSlow(:,2,:))];
    medCorrInfraGamma = [medCorrInfraGamma; squeeze(mVars.medCorrInfraSlow(:,3,:))];
    medCorrInfraWB    = [medCorrInfraWB; squeeze(mVars.medCorrInfraSlow(:,4,:))];

    % Dividing pairwise correlations based on superficial and deeper layers
    % of cortex
    superMeanPairCorrAlpha = [superMeanPairCorrAlpha; squeeze(mVars.superMeanCorr(:,1,:))];
    superMeanPairCorrBeta  = [superMeanPairCorrBeta; squeeze(mVars.superMeanCorr(:,2,:))];
    superMeanPairCorrGamma = [superMeanPairCorrGamma; squeeze(mVars.superMeanCorr(:,3,:))];
    superMeanPairCorrWB    = [superMeanPairCorrWB; squeeze(mVars.superMeanCorr(:,4,:))];

    superMedPairCorrAlpha = [superMedPairCorrAlpha; squeeze(mVars.superMedCorr(:,1,:))];
    superMedPairCorrBeta  = [superMedPairCorrBeta; squeeze(mVars.superMedCorr(:,2,:))];
    superMedPairCorrGamma = [superMedPairCorrGamma; squeeze(mVars.superMedCorr(:,3,:))];
    superMedPairCorrWB    = [superMedPairCorrWB; squeeze(mVars.superMedCorr(:,4,:))];

    superMeanCorrEnvAlpha = [superMeanCorrEnvAlpha; squeeze(mVars.superMeanEnvelopeCorr(:,1,:))];
    superMeanCorrEnvBeta  = [superMeanCorrEnvBeta; squeeze(mVars.superMeanEnvelopeCorr(:,2,:))];
    superMeanCorrEnvGamma = [superMeanCorrEnvGamma; squeeze(mVars.superMeanEnvelopeCorr(:,3,:))];
    superMeanCorrEnvWB    = [superMeanCorrEnvWB; squeeze(mVars.superMeanEnvelopeCorr(:,4,:))];

    superMedCorrEnvAlpha = [superMedCorrEnvAlpha; squeeze(mVars.superMedEnvelopeCorr(:,1,:))];
    superMedCorrEnvBeta  = [superMedCorrEnvBeta; squeeze(mVars.superMedEnvelopeCorr(:,2,:))];
    superMedCorrEnvGamma = [superMedCorrEnvGamma; squeeze(mVars.superMedEnvelopeCorr(:,3,:))];
    superMedCorrEnvWB    = [superMedCorrEnvWB; squeeze(mVars.superMedEnvelopeCorr(:,4,:))];

    superMeanCorrInfraAlpha = [superMeanCorrInfraAlpha; squeeze(mVars.superMeanInfraSlowCorr(:,1,:))];
    superMeanCorrInfraBeta  = [superMeanCorrInfraBeta; squeeze(mVars.superMeanInfraSlowCorr(:,2,:))];
    superMeanCorrInfraGamma = [superMeanCorrInfraGamma; squeeze(mVars.superMeanInfraSlowCorr(:,3,:))];
    superMeanCorrInfraWB    = [superMeanCorrInfraWB; squeeze(mVars.superMeanInfraSlowCorr(:,4,:))];

    superMedCorrInfraAlpha = [superMedCorrInfraAlpha; squeeze(mVars.superMedInfraSlowCorr(:,1,:))];
    superMedCorrInfraBeta  = [superMedCorrInfraBeta; squeeze(mVars.superMedInfraSlowCorr(:,2,:))];
    superMedCorrInfraGamma = [superMedCorrInfraGamma; squeeze(mVars.superMedInfraSlowCorr(:,3,:))];
    superMedCorrInfraWB    = [superMedCorrInfraWB; squeeze(mVars.superMedInfraSlowCorr(:,4,:))];

    deepMeanPairCorrAlpha = [deepMeanPairCorrAlpha; squeeze(mVars.deepMeanCorr(:,1,:))];
    deepMeanPairCorrBeta  = [deepMeanPairCorrBeta; squeeze(mVars.deepMeanCorr(:,2,:))];
    deepMeanPairCorrGamma = [deepMeanPairCorrGamma; squeeze(mVars.deepMeanCorr(:,3,:))];
    deepMeanPairCorrWB    = [deepMeanPairCorrWB; squeeze(mVars.deepMeanCorr(:,4,:))];

    deepMedPairCorrAlpha = [deepMedPairCorrAlpha; squeeze(mVars.deepMedCorr(:,1,:))];
    deepMedPairCorrBeta  = [deepMedPairCorrBeta; squeeze(mVars.deepMedCorr(:,2,:))];
    deepMedPairCorrGamma = [deepMedPairCorrGamma; squeeze(mVars.deepMedCorr(:,3,:))];
    deepMedPairCorrWB    = [deepMedPairCorrWB; squeeze(mVars.deepMedCorr(:,4,:))];

    deepMeanCorrEnvAlpha = [deepMeanCorrEnvAlpha; squeeze(mVars.deepMeanEnvelopeCorr(:,1,:))];
    deepMeanCorrEnvBeta  = [deepMeanCorrEnvBeta; squeeze(mVars.deepMeanEnvelopeCorr(:,2,:))];
    deepMeanCorrEnvGamma = [deepMeanCorrEnvGamma; squeeze(mVars.deepMeanEnvelopeCorr(:,3,:))];
    deepMeanCorrEnvWB    = [deepMeanCorrEnvWB; squeeze(mVars.deepMeanEnvelopeCorr(:,4,:))];

    deepMedCorrEnvAlpha = [deepMedCorrEnvAlpha; squeeze(mVars.deepMedEnvelopeCorr(:,1,:))];
    deepMedCorrEnvBeta  = [deepMedCorrEnvBeta; squeeze(mVars.deepMedEnvelopeCorr(:,2,:))];
    deepMedCorrEnvGamma = [deepMedCorrEnvGamma; squeeze(mVars.deepMedEnvelopeCorr(:,3,:))];
    deepMedCorrEnvWB    = [deepMedCorrEnvWB; squeeze(mVars.deepMedEnvelopeCorr(:,4,:))];

    deepMeanCorrInfraAlpha = [deepMeanCorrInfraAlpha; squeeze(mVars.deepMeanInfraSlowCorr(:,1,:))];
    deepMeanCorrInfraBeta  = [deepMeanCorrInfraBeta; squeeze(mVars.deepMeanInfraSlowCorr(:,2,:))];
    deepMeanCorrInfraGamma = [deepMeanCorrInfraGamma; squeeze(mVars.deepMeanInfraSlowCorr(:,3,:))];
    deepMeanCorrInfraWB    = [deepMeanCorrInfraWB; squeeze(mVars.deepMeanInfraSlowCorr(:,4,:))];

    deepMedCorrInfraAlpha = [deepMedCorrInfraAlpha; squeeze(mVars.deepMedInfraSlowCorr(:,1,:))];
    deepMedCorrInfraBeta  = [deepMedCorrInfraBeta; squeeze(mVars.deepMedInfraSlowCorr(:,2,:))];
    deepMedCorrInfraGamma = [deepMedCorrInfraGamma; squeeze(mVars.deepMedInfraSlowCorr(:,3,:))];
    deepMedCorrInfraWB    = [deepMedCorrInfraWB; squeeze(mVars.deepMedInfraSlowCorr(:,4,:))];

    for iList = 1:size(mVars.marginalA_B_temp,2)
        if ~isempty(mVars.marginalA_B_temp{iList})
            % Normalize the vector
            numRows = min([size(mVars.marginalA_B_temp{iList},1) size(mVars.marginalB_A_temp{iList},1)]);
            if numRows>=20
                numRows = 20;
            else
                numRows = min([size(mVars.marginalA_B_temp{iList},1) size(mVars.marginalB_A_temp{iList},1)]);
            end
            alphaVarsAB(:,:,rowIdx) = squeeze([mVars.marginalA_B_temp{iList}(1:numRows,1,:)./max(mVars.marginalA_B_temp{iList}(1:numRows,1,:)) ;zeros(20-numRows,1,2)]);
            alphaVarsBA(:,:,rowIdx) = squeeze([mVars.marginalB_A_temp{iList}(1:numRows,1,:)./max(mVars.marginalB_A_temp{iList}(1:numRows,1,:)) ;zeros(20-numRows,1,2)]);

            gammaVarsAB(:,:,rowIdx) = squeeze([mVars.marginalA_B_temp{iList}(1:numRows,3,:)./max(mVars.marginalA_B_temp{iList}(1:numRows,3,:)) ;zeros(20-numRows,1,2)]);
            gammaVarsBA(:,:,rowIdx) = squeeze([mVars.marginalB_A_temp{iList}(1:numRows,3,:)./max(mVars.marginalB_A_temp{iList}(1:numRows,3,:)) ;zeros(20-numRows,1,2)]);

            one_oneMapAlpha(:,:,rowIdx) = squeeze([mVars.one_oneCorr_temp{iList}(1:numRows,1,:)./max(mVars.one_oneCorr_temp{iList}(1:numRows,1,:)) ;zeros(20-numRows,1,2)]);
            one_oneMapGamma(:,:,rowIdx) = squeeze([mVars.one_oneCorr_temp{iList}(1:numRows,3,:)./max(mVars.one_oneCorr_temp{iList}(1:numRows,3,:)) ;zeros(20-numRows,1,2)]);
            rowIdx = rowIdx+1;
        else
            removeIdx = [removeIdx rowIdx];
            rowIdx = rowIdx+1;

        end
    end
end
badDatIdx = [badDatIdx; 77;78;79]; % Remove anesthesia controls from data

%% Plot pairwise correlations with FC and distance
cVal = {'b';'r';[0 0.4470 0.7410];'m';'k'};

for iBand = 3%1:4
    fcTemp = fcVals; distTemp = dist; pairTypeTemp = pairType; patchTemp = patchInfo;
    clear meanPairCorr meanCorrEnv meanCorrInfra superPairCorr superEnvCorr superInfraCorr deepPairCorr deepEnvCorr deepInfraCorr bandName
    switch iBand
        case 1
            bandName = 'Alpha';
            meanPairCorr = meanPairCorrAlpha;
            meanCorrEnv  = meanCorrEnvAlpha;
            meanCorrInfra = meanCorrInfraAlpha;

            superPairCorr = superMeanPairCorrAlpha;
            superEnvCorr = superMeanCorrEnvAlpha;
            superInfraCorr = superMeanCorrInfraAlpha;

            deepPairCorr = deepMeanPairCorrAlpha;
            deepEnvCorr = deepMeanCorrEnvAlpha;
            deepInfraCorr = deepMeanCorrInfraAlpha;

        case 2
            bandName = 'Beta';
            meanPairCorr = meanPairCorrBeta;
            meanCorrEnv  = meanCorrEnvBeta;
            meanCorrInfra = meanCorrInfraBeta;

            superPairCorr = superMeanPairCorrBeta;
            superEnvCorr = superMeanCorrEnvBeta;
            superInfraCorr = superMeanCorrInfraBeta;

            deepPairCorr = deepMeanPairCorrBeta;
            deepEnvCorr = deepMeanCorrEnvBeta;
            deepInfraCorr = deepMeanCorrInfraBeta;

        case 3
            bandName = 'Gamma';
            meanPairCorr = meanPairCorrGamma;
            meanCorrEnv  = meanCorrEnvGamma;
            meanCorrInfra = meanCorrInfraGamma;

            superPairCorr = superMeanPairCorrGamma;
            superEnvCorr = superMeanCorrEnvGamma;
            superInfraCorr = superMeanCorrInfraGamma;

            deepPairCorr = deepMeanPairCorrGamma;
            deepEnvCorr = deepMeanCorrEnvGamma;
            deepInfraCorr = deepMeanCorrInfraGamma;

        case 4
            bandName = 'Wideband';
            meanPairCorr = meanPairCorrWB;
            meanCorrEnv  = meanCorrEnvWB;
            meanCorrInfra = meanCorrInfraWB;

            superPairCorr = superMeanPairCorrWB;
            superEnvCorr = superMeanCorrEnvWB;
            superInfraCorr = superMeanCorrInfraWB;

            deepPairCorr = deepMeanPairCorrWB;
            deepEnvCorr = deepMeanCorrEnvWB;
            deepInfraCorr = deepMeanCorrInfraWB;

    end

    meanPairCorr(badDatIdx,:) = []; meanCorrEnv(badDatIdx,:) = []; meanCorrInfra(badDatIdx,:) = [];
    superPairCorr(badDatIdx,:) = []; superEnvCorr(badDatIdx,:) = []; superInfraCorr(badDatIdx,:) = [];
    deepPairCorr(badDatIdx,:) = []; deepEnvCorr(badDatIdx,:) = []; deepInfraCorr(badDatIdx,:) = [];

    fcTemp(badDatIdx) = []; distTemp(badDatIdx)= []; pairTypeTemp(badDatIdx) = []; patchTemp(badDatIdx) = []; patchTemp = logical(patchTemp);
    ssPairs = find(strcmp(pairTypeTemp,'SS')); smPairs = find(strcmp(pairTypeTemp,'SM')); mmPairs = find(strcmp(pairTypeTemp,'MM'));

    for iRef = 1%:2
        if iRef== 1; refTitle = 'No referencing'; else; refTitle = 'Avg Reference'; end

        plotPairCorrs(fcTemp,meanPairCorr(:,iRef),meanCorrEnv(:,iRef),meanCorrInfra(:,iRef),...
            superPairCorr(:,iRef), superEnvCorr(:,iRef),superInfraCorr(:,iRef),...
            deepPairCorr(:,iRef),deepEnvCorr(:,iRef),deepInfraCorr(:,iRef),...
            'Functional Connectivity',cVal{iBand},'poly1');
        sgtitle(['Pairwise correlations vs Functional Connectivity for ' bandName ' for ' refTitle]);

        % Pairwise vs distance
        plotPairCorrs(distTemp,meanPairCorr(:,iRef),meanCorrEnv(:,iRef),meanCorrInfra(:,iRef),...
            superPairCorr(:,iRef), superEnvCorr(:,iRef),superInfraCorr(:,iRef),...
            deepPairCorr(:,iRef),deepEnvCorr(:,iRef),deepInfraCorr(:,iRef),...
            'Distance ',cVal{iBand},'poly2');
        sgtitle(['Pairwise correlations vs Distance for ' bandName ' for ' refTitle]);

        % Pairwise vs Functional connectivity for SS, SM, MM pairs
        plotPairCorrs(fcTemp(ssPairs),meanPairCorr(ssPairs,iRef),meanCorrEnv(ssPairs,iRef),meanCorrInfra(ssPairs,iRef),...
            superPairCorr(ssPairs,iRef), superEnvCorr(ssPairs,iRef),superInfraCorr(ssPairs,iRef),...
            deepPairCorr(ssPairs,iRef),deepEnvCorr(ssPairs,iRef),deepInfraCorr(ssPairs,iRef),...
            'Functional Connectivity',cVal{iBand},'poly1'); % Sensory-Sensory pairs
        sgtitle(['Pairwise correlations vs Functional Connectivity for Sensory-sensory pairs for' bandName ' for ' refTitle]);

        plotPairCorrs(fcTemp(smPairs),meanPairCorr(smPairs,iRef),meanCorrEnv(smPairs,iRef),meanCorrInfra(smPairs,iRef),...
            superPairCorr(smPairs,iRef), superEnvCorr(smPairs,iRef),superInfraCorr(smPairs,iRef),...
            deepPairCorr(smPairs,iRef),deepEnvCorr(smPairs,iRef),deepInfraCorr(smPairs,iRef),...
            'Functional Connectivity',cVal{iBand},'poly1'); % Sensory-Motor pairs
        sgtitle(['Pairwise correlations vs Functional Connectivity for Sensory-motor pairs for  ' bandName ' for ' refTitle]);

        plotPairCorrs(fcTemp(mmPairs),meanPairCorr(mmPairs,iRef),meanCorrEnv(mmPairs,iRef),meanCorrInfra(mmPairs,iRef),...
            superPairCorr(mmPairs,iRef), superEnvCorr(mmPairs,iRef),superInfraCorr(mmPairs,iRef),...
            deepPairCorr(mmPairs,iRef),deepEnvCorr(mmPairs,iRef),deepInfraCorr(mmPairs,iRef),...
            'Functional Connectivity',cVal{iBand},'poly1'); % Motor-Motor pairs
        sgtitle(['Pairwise correlations vs Functional Connectivity for motor-motor pairs for ' bandName ' for ' refTitle]);

        % Pairwise vs Functional connectivity for sites in patch and out
        % of patch
        plotPairCorrs(fcTemp(patchTemp),meanPairCorr(patchTemp,iRef),meanCorrEnv(patchTemp,iRef),meanCorrInfra(patchTemp,iRef),...
            superPairCorr(patchTemp,iRef), superEnvCorr(patchTemp,iRef),superInfraCorr(patchTemp,iRef),...
            deepPairCorr(patchTemp,iRef),deepEnvCorr(patchTemp,iRef),deepInfraCorr(patchTemp,iRef),...
            'Functional Connectivity',cVal{iBand},'poly1');
        sgtitle(['Pairwise correlations vs Functional Connectivity for site pairs in patch for ' bandName ' for ' refTitle]);


        plotPairCorrs(fcTemp(~patchTemp),meanPairCorr(~patchTemp,iRef),meanCorrEnv(~patchTemp,iRef),meanCorrInfra(~patchTemp,iRef),...
            superPairCorr(~patchTemp,iRef), superEnvCorr(~patchTemp,iRef),superInfraCorr(~patchTemp,iRef),...
            deepPairCorr(~patchTemp,iRef),deepEnvCorr(~patchTemp,iRef),deepInfraCorr(~patchTemp,iRef),...
            'Functional Connectivity',cVal{iBand},'poly1');
        sgtitle(['Pairwise correlations vs Functional Connectivity for site pairs out of patch for ' bandName ' for ' refTitle]);
    end
end

%% Normalize and plot the infraslow pairwise correlations inside and outside the patch
fcTemp = fcVals; distTemp = dist; pairTypeTemp = pairType; patchTemp = patchInfo;
fcTemp(badDatIdx) = []; distTemp(badDatIdx)= []; pairTypeTemp(badDatIdx) = []; patchTemp(badDatIdx) = []; patchTemp = logical(patchTemp);

one_oneMapAlphaTemp = one_oneMapAlpha; one_oneMapGammaTemp = one_oneMapGamma;
one_oneMapAlphaTemp(:,:,badDatIdx) = []; one_oneMapGammaTemp(:,:,badDatIdx)= [];


edgeInfoFC = []; edgeInfoFCP = []; edgeInfoFCNP= [];   edgeInfoDist = []; edgeInfoDistP = []; edgeInfoDistNP = [];

for iRef = 1%:2
    switch iRef
        case 1
            refName = 'No referencing';
        case 2
            refName = 'Avg referencing';
    end

    for iType =2% 1:2
        clear sortVar sortInd patchNew alphaTemp gammaTemp one_oneMapAlphaPTemp one_oneMapAlphaNPTemp one_oneMapAlphaNPTemp one_oneMapGammaNPTemp...
            label edges medColValsAlpha medColValsGamma medColValsAlphaNoPatch medColValsGammaNoPatch medColValsAlphaPatch medColValsGammaPatch
        
        switch iType
            case 1 % FC
                [sortVar, sortInd] = sort(fcTemp);
                patchNew           = patchTemp(sortInd);
                edges              = -0.4:0.2:0.8;                
                label              = 'FC ';

            case 2
                [sortVar,sortInd] = sort(distTemp);
                patchNew          = patchTemp(sortInd);
                edges             = 0:2:18;                
                label             = 'Distance ';
        end

        edgeInfoVar = []; edgeInfoP = [];  edgeInfoNP = [];

        alphaTemp = squeeze(one_oneMapAlphaTemp(:,iRef,sortInd));
        gammaTemp = squeeze(one_oneMapGammaTemp(:,iRef,sortInd));

        one_oneMapAlphaPTemp  = alphaTemp(:,patchNew);
        one_oneMapAlphaNPTemp = alphaTemp(:,~patchNew);

        one_oneMapGammaPTemp  = gammaTemp(:,patchNew);
        one_oneMapGammaNPTemp = gammaTemp(:,~patchNew);

        % Sorting
        [y,edgeDist] = discretize(sortVar,edges);
        binVals      = unique(y);

        for iBin = 1:length(edges)
            colVals = find(y == iBin);

            if isempty(colVals)
                medColValsAlpha(:,iBin) = NaN(20,1);
                medColValsGamma(:,iBin) = NaN(20,1);
                continue;
            end

            medColValsAlpha(:,iBin) = median(alphaTemp(:,colVals),2,'omitnan');
            medColValsGamma(:,iBin) = median(gammaTemp(:,colVals),2,'omitnan');
            edgeInfoVar = [edgeInfoVar edgeDist(iBin)];
        end

        % Sorting - within a patch
        clear edgeDist
        [yPatch,edgeDist] = discretize(sortVar(patchNew),edges);
        binPatch = unique(yPatch);

        for iBin = 1:length(edges)
            colVals = find(yPatch == iBin);
            if isempty(colVals)
                medColValsAlphaPatch(:,iBin) = NaN(20,1);
                medColValsGammaPatch(:,iBin) = NaN(20,1);
                continue;
            end
            medColValsAlphaPatch(:,iBin) = median(one_oneMapAlphaPTemp(:,colVals),2,'omitnan');
            medColValsGammaPatch(:,iBin) = median(one_oneMapGammaPTemp(:,colVals),2,'omitnan');
            edgeInfoP = [edgeInfoP edgeDist(iBin)];
        end

        % Sorting with FC - outside a patch
        clear edgeDist
        [yNoPatch, edgeDist] = discretize(sortVar(~patchNew),edges);
        binNoPatch = unique(yNoPatch);

        for iBin = 1:length(edges)%max(binNoPatch)
            colVals = find(yNoPatch == iBin);
            if isempty(colVals)
                medColValsAlphaFCNoPatch(:,iBin) = NaN(20,1);
                medColValsGammaFCNoPatch(:,iBin) = NaN(20,1);
                continue;
            end
            medColValsAlphaNoPatch(:,iBin) = median(one_oneMapAlphaNPTemp(:,colVals),2,'omitnan');
            medColValsGammaNoPatch(:,iBin) = median(one_oneMapGammaNPTemp(:,colVals),2,'omitnan');
            edgeInfoNP = [edgeInfoNP edgeDist(iBin)];
        end

        medColValsAlpha(:,all(isnan(medColValsAlpha),1))=[];
        if find(all(medColValsAlpha==0)); removeCols = find(all(medColValsAlpha==0)); medColValsAlpha(:,removeCols) =[]; end 

        medColValsGamma(:,all(isnan(medColValsGamma),1))=[];
        if find(all(medColValsGamma==0)); removeCols = find(all(medColValsGamma==0)); medColValsGamma(:,removeCols) =[]; end 

        medColValsAlphaNoPatch(:,all(isnan(medColValsAlphaNoPatch),1))= [];
        if find(all(medColValsAlphaNoPatch==0)); removeCols = find(all(medColValsAlphaNoPatch==0)); medColValsAlphaNoPatch(:,removeCols) =[]; end 

        medColValsGammaNoPatch(:,all(isnan(medColValsGammaNoPatch),1))=[];
        if find(all(medColValsGammaNoPatch==0)); removeCols = find(all(medColValsGammaNoPatch==0)); medColValsGammaNoPatch(:,removeCols) =[]; end 

        medColValsAlphaPatch(:,all(isnan(medColValsAlphaPatch),1)) = [];
        if find(all(medColValsAlphaPatch==0)); removeCols = find(all(medColValsAlphaPatch==0)); medColValsAlphaPatch(:,removeCols) =[]; end 

        medColValsGammaPatch(:,all(isnan(medColValsGammaPatch),1)) = [];
        if find(all(medColValsGammaPatch==0)); removeCols = find(all(medColValsGammaPatch==0)); medColValsGammaPatch(:,removeCols) =[]; end      


        % Plot all sorted data for alpha
        figure;subplot(131); imagesc(imgaussfilt(movmean(medColValsAlpha,2,1),1));
        clim([0 1]); colorbar; colormap jet; shading interp; axis square;
        xticks(1:size(medColValsAlpha,2));
        xticklabels(floor(edgeInfoVar.*100)./100); title(label);

        subplot(132); imagesc(imgaussfilt(movmean(medColValsAlphaPatch,2,1),1));
        clim([0 1]); colorbar; colormap jet;shading interp;axis square;
        xticks(1:size(medColValsAlphaPatch,2));
        xticklabels(floor(edgeInfoP.*100)./100);
        title([label ' In patch ' ]);

        subplot(133); imagesc(imgaussfilt(movmean(medColValsAlphaNoPatch,2,1),1));
        clim([0 1]); colorbar; colormap jet;shading interp;axis square;
        xticks(1:size(medColValsAlphaNoPatch,2));
        xticklabels(floor(edgeInfoNP.*100)./100);
        title([label ' Out of patch ' ]);
        sgtitle(['Alpha band ' refName]);

        % Plot all sorted data for gamma
        figure;subplot(131); imagesc(imgaussfilt(movmean(medColValsGamma,2,1),1));
        clim([0.3 1]); colorbar; colormap jet;shading interp;axis square;
        xticks(1:size(medColValsGamma,2));
        xticklabels(floor(edgeInfoVar.*100)./100);
        title( label );

        subplot(132); imagesc(imgaussfilt(movmean(medColValsGammaPatch,2,1),1));
        clim([0.3 1]); colorbar; colormap jet;shading interp;
        xticks(1:size(medColValsGammaPatch,2));axis square;
        xticklabels(floor(edgeInfoP.*100)./100);
        title([label '  In patch ']);

        subplot(133); imagesc(imgaussfilt(movmean(medColValsGammaNoPatch,2,1),1));
        clim([0.3 1]); colorbar; colormap jet;shading interp;axis square;
        xticks(1:size(medColValsGammaNoPatch,2));
        xticklabels(floor(edgeInfoNP.*100)./100);
        title([label 'Out of patch ']);
        sgtitle(['Gamma band ' refName]);
    end

end
%% Models explaining pairwise corr and FC
infraCorrTemp =  meanCorrInfraGamma(:,1); infraCorrTemp(badDatIdx)=[];
X4 = [ones(size(distTemp))  infraCorrTemp distTemp ]; 
[b,bint,r,rint,stats] = regress(fcTemp,X4);

%%  Function to plot all scatter plots
function plotPairCorrs(x,meanPairCorr,meanCorrEnv,meanCorrInfra,superPairCorr, superEnvCorr,superInfraCorr,...
    deepPairCorr,deepEnvCorr,deepInfraCorr,xLabel,cVal,fitType)

figure;
subplot(331);scatter(x,meanPairCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,meanPairCorr,fitType,cVal); xlabel(xLabel); ylabel('LFP Time series correlation'); title('All Channels');

subplot(332);scatter(x,meanCorrEnv,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,meanCorrEnv,fitType,cVal);xlabel(xLabel); ylabel('Instantaneous power correlation');title('All Channels');

subplot(333);scatter(x,meanCorrInfra,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,meanCorrInfra,fitType,cVal); xlabel(xLabel); ylabel('Infraslow power correlation');title('All Channels');

subplot(334);scatter(x,superPairCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,superPairCorr,fitType,cVal);xlabel(xLabel); ylabel('LFP Time series correlation');title('Superficial Channels');

subplot(335);scatter(x,superEnvCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,superEnvCorr,fitType,cVal);xlabel(xLabel); ylabel('Instantaneous power correlation');title('Superficial Channels');

subplot(336);scatter(x,superInfraCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,superInfraCorr,fitType,cVal);xlabel(xLabel); ylabel('Infraslow power correlation');title('Superficial Channels');

subplot(337);scatter(x,deepPairCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,deepPairCorr,fitType,cVal);xlabel(xLabel); ylabel('LFP Time series correlation'); title('Deep Channels');

subplot(338);scatter(x,deepEnvCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,deepEnvCorr,fitType,cVal);xlabel(xLabel); ylabel('Instantaneous power correlation');title('Deep Channels');

subplot(339);scatter(x,deepInfraCorr,35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal);hold on;
plotFit(x,deepInfraCorr,fitType,cVal); xlabel(xLabel); ylabel('Infraslow power correlation');title('Deep Channels');
end
%% Function to fit a line and plot mean+Std dev as patch
function plotFit(x,var,fitType,colorVals)
include = ~isnan(var);
xFit     = linspace(min(x(include)),max(x(include)),1000);
[f,gof]  = fit(x(include),var(include),fitType);
c        = coeffvalues(f);
r2       = gof.rsquare;

if strcmp(fitType,'poly1')
    fitLine  = c(1)*(xFit) +c(2);
elseif strcmp(fitType,'poly2')
    fitLine = c(1)*(xFit.^2) + c(2)*(xFit) +c(3);
end

plot(xFit,fitLine,'Color',colorVals,'LineStyle','--','LineWidth',2);
text(-0.2,0.8,['R^2: ' num2str(r2)],'Color',colorVals);
patch([xFit fliplr(xFit)],[fitLine-2*std(fitLine) fliplr(fitLine +2*std(fitLine))],colorVals,'FaceAlpha',0.2,'EdgeColor','none');
ylim([-0.5 1]);
if max(x)<=1; xlim([-0.5 1]); else; xlim([0 15]); end
end