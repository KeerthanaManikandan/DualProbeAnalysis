%% ephysDualProbeRS_v7
% This function performs analysis on LFP recorded simultaneously from two linear electrode arrays.
% October 15, 2025 - Keerthana Manikandan
% This code performs the following for ONE MONKEY (for combined analyis,
% check this script- combinedAnalysisEphysDualProbeRS.m)
% Check ephysDualProbeRS_v6 for previous iterations

% Set paths
clear;clc
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\DualProbe']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\spectral_analysis\continuous\dupes']));
clc;

%% Initialize all relevant variables
clear allDates dateFileNumAll serverPath refDate refDir refImageName datFileNameAll chInCortexProbeA ...
    chInCortexProbeB probeLabelA probeLabelB anesthesiaLevels heartRate patchInfo pairClass

fs        = 1e3; % Sampling frequency
thetaBand = [6 8];    [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters
alphaBand = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass');
betaBand  = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass');

hemisphere  = 'Left';
saveFigureFlag = 0;


chOutCortex    = 1:3;
chDeep         = 30:32;

iM = 1; % 1 - Charlie Sheen, 2 - Whiskey
switch iM
    case 1
        monkeyName = 'CharlieSheen';
        goodRuns = [1 0 1 1 1 1 1 1 0 1 1 NaN NaN NaN NaN NaN; ...
            NaN NaN NaN NaN 1 0 1 1 1 1 NaN NaN NaN NaN NaN NaN; ...
            NaN NaN NaN NaN 1 1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN; ...
            1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN NaN NaN; ...
            1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';

    case 2
        monkeyName = 'Whiskey';
        goodRuns = [1 1 1 1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; ...
            1 1 1 0 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
            1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN;...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 NaN NaN]';
end

% Load all necessary variables
[allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate,patchInfo,pairClass] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']);

% Get the connectivity and distance between pairs
clear distSites connSites greenMapRef
[distSites,connVals,refSites,movSites,greenMapRef] = ...
    getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);


distSitesAll = NaN(size(goodRuns)); connValsAll = NaN(size(goodRuns));
heartRateValsAll = NaN(size(goodRuns)); anesthesiaValsAll = NaN(size(goodRuns));


for iDate =1:size(allDates,1)
    clear datFileNum
    datFileNum        = datFileNumAll{iDate,1};
    distSitesAll(datFileNum,iDate) = distSites{iDate,1}(datFileNum);
    connValsAll(datFileNum,iDate)  = squeeze(mean(connVals{iDate,1}(:,:,datFileNum),[1,2],'omitnan'));

    anesthesiaValsAll(datFileNum,iDate) = anesthesiaLevels{iDate,1}(datFileNum);
    heartRateValsAll(datFileNum,iDate)  = heartRate{iDate,1}(datFileNum);
end

disp(['Obtained/retrieved distance between probes, connectivity values, heart rate, anesthesia levels for ' monkeyName]);

% Get all dual probe data
[allProbeData,allBadTimes,badElecA,badElecB,estChInCortexA,estChInCortexB] = getAllDualProbeData(monkeyName,hemisphere,allDates, ...
    datFileNameAll, datFileNumAll,serverPath,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag );

%% Get all pairwise correlations
bandLabels = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking'};
timeLabels = {'Time series','Power','Infraslow'};
[allVars] = getDualProbeCorrelations(monkeyName, hemisphere, allDates, datFileNumAll,allProbeData,allBadTimes,...
    badElecA,badElecB,estChInCortexA,estChInCortexB,connValsAll,distSitesAll,goodRuns);

% Check plotting...


%% Perform Phase amplitude coupling
% Get the PAC comodulogram
% Set the amplitude range
clear modIdxAllA2B modIdxAllB2A modIdxAllA2A modIdxAllB2B
gammaRange = gammaBand(1):5:gammaBand(2);
gammaRange(gammaRange>=60 & gammaRange<=65)=[];

% Get the phase range
lowFreqRange = 6:2:30;

nHigh    = size(gammaRange,2);
nLow     = size(lowFreqRange,2);

% Initialize window and step size
winSize  = 2e3;
stepSize = 2e3;

% Get the phase amplitude coupling for the recordings...
for iDate = 2:5%size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        fileMat = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_2sec_' num2str(fileNum) '.mat'];
        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);

        if chA(1)== 0 || chB(1)==0; continue; end

        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_2sec_' num2str(fileNum) '.mat'],'file') %(~ismember('modIdxAllA2BCircleT', who('-file', fileMat))) || 1%
            clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);


            [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
                allProbeData{fileNum,iDate}.probe1Ch,allProbeData{fileNum,iDate}.probe2Ch,badElecA{fileNum,iDate},...
                badElecB{fileNum,iDate},allBadTimes{fileNum,iDate},estChInCortexA{iDate}(iRun,:),...
                estChInCortexB{iDate}(iRun,:),gammaRange,lowFreqRange);

            % Get comodulogram
            nChan    = min(size(amplitudeA,3),size(amplitudeB,3));
            dataLen  = size(amplitudeA,2);


            % Calculate windows
            winStart = 1:stepSize:(dataLen-winSize);
            winEnd   = winStart+winSize-1;
            nWin     = numel(winStart);

            % Preallocating modulation indices
            modA2B = NaN(nWin,nHigh,nLow,nChan,'single');
            modB2A = NaN(nWin,nHigh,nLow,nChan,'single');
            modA2A = NaN(nWin,nHigh,nLow,nChan,'single');
            modB2B = NaN(nWin,nHigh,nLow,nChan,'single');

            % Avoid overhead
            lowFreqAconst  = parallel.pool.Constant(phaseA);
            lowFreqBconst  = parallel.pool.Constant(phaseB);
            highFreqAconst = parallel.pool.Constant(amplitudeA);
            highFreqBconst = parallel.pool.Constant(amplitudeB);


            tic;
            parfor iHigh = 1:nHigh
                % Local copies
                lowA = lowFreqAconst.Value;
                lowB = lowFreqBconst.Value;
                highA = highFreqAconst.Value;
                highB = highFreqBconst.Value;

                modA2BT = NaN(nWin,nLow,nChan,'single');
                modB2AT = NaN(nWin,nLow,nChan,'single');
                modA2AT = NaN(nWin,nLow,nChan,'single');
                modB2BT = NaN(nWin,nLow,nChan,'single');

                for iWin = 1:nWin
                    idx = winStart(iWin):winEnd(iWin);
                    for iLow = 1:nLow
                        [modA2BT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)));
                        [modB2AT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)));
                        [modA2AT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)));
                        [modB2BT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)));
                    end
                end

                modA2B(:,iHigh,:,:) = modA2BT;
                modB2A(:,iHigh,:,:) = modB2AT;
                modA2A(:,iHigh,:,:) = modA2AT;
                modB2B(:,iHigh,:,:) = modB2BT;

            end
            toc;

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_2sec_' num2str(fileNum) '.mat'],...
                'modA2B','modB2A','modA2A','modB2B');
            clear lowFreqAconst lowFreqBconst highFreqAconst highFreqBconst

            modIdxAllA2B2sec{iRun,iDate} = modA2B;
            modIdxAllB2A2sec{iRun,iDate} = modB2A;
            modIdxAllA2A2sec{iRun,iDate} = modA2A;
            modIdxAllB2B2sec{iRun,iDate} = modB2B;

        else
            clear vars;
            vars = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_2sec_' num2str(fileNum) '.mat']);
            modIdxAllA2B2sec{iRun,iDate} = vars.modA2B;
            modIdxAllB2A2sec{iRun,iDate} = vars.modB2A;
            modIdxAllA2A2sec{iRun,iDate} = vars.modA2A;
            modIdxAllB2B2sec{iRun,iDate} = vars.modB2B;

        end
    end
end

%% 
a2b2Sec = modA2B;
b2a2Sec = modB2A;
a2a2Sec = modA2A;
b2b2Sec = modB2B;

a2b5Sec = modA2B;
b2a5Sec = modB2A;
a2a5Sec = modA2A;
b2b5Sec = modB2B;

a2b10Sec = modA2B;
b2a10Sec = modB2A;
a2a10Sec = modA2A;
b2b10Sec = modB2B;

figure; 
subplot(441); imagesc(lowFreqRange,gammaRange,squeeze(median(modIdxAllB2B{iRun,iDate},[1 4],'omitnan'))); 
axis square; colorbar
title('100 s window');
set(gca,'YDir','normal');



%% Plot the comodulogram - edit this
% Average across all channels for the 4 combinations...
avgModA2B = squeeze(mean(modIdxAllA2B{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
avgModB2A = squeeze(mean(modIdxAllB2A{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
avgModA2A = squeeze(mean(modIdxAllA2A{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
avgModB2B = squeeze(mean(modIdxAllB2B{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));

figure;
subplot(221); contourf(lowFreqRange,gammaRange,median(avgModA2B,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet; title('A--B');
subplot(222); contourf(lowFreqRange,gammaRange,median(avgModB2A,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;title('B--A');
subplot(223); contourf(lowFreqRange,gammaRange,median(avgModA2A,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;title('A--A');
subplot(224); contourf(lowFreqRange,gammaRange,median(avgModB2B,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;title('B--B');

% Plot it for all channels
for iPlot = 1:4
    switch iPlot
        case 1
            plotVar = avgModA2B;
            figTitle = 'A--B';
        case 2
            plotVar = avgModB2A;
            figTitle = 'B--A';
        case 3
            plotVar = avgModA2A;
            figTitle = 'A--A';
        case 4
            plotVar = avgModB2B;
            figTitle = 'B--B';
    end

    figure;
    for iCh = 1: size(plotVar,3)
        subplot(size(plotVar,3),1,iCh)
        contourf(lowFreqRange,gammaRange,plotVar(:,:,iCh),'lines','none'); colormap jet; clim([0 1e-4])
        if iCh~=size(plotVar,3)
            xticklabels({}); yticklabels({});
        end
    end
    sgtitle(figTitle);
end

% Plot for compartments....
%% Circular shifting the comodulogram


%% Shuffling the windows in increments
% winSize = [10 20 30 40 50 100].*1e3;
nHigh = size(gammaRange,2);
nLow  = size(lowFreqRange,2);


% shiftLen = [1 5 10 20 50 100].*1e3;
% nShift   = length(shiftLen);

for iDate = 2:5
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);
        if chA(1)== 0 || chB(1)==0; continue; end
        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate...
                '\Electrophysiology\modulogramCtrl_' num2str(fileNum) '.mat'],'file')
            clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

            [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
                allProbeData{fileNum,iDate}.probe1Ch,allProbeData{fileNum,iDate}.probe2Ch,badElecA{fileNum,iDate},...
                badElecB{fileNum,iDate},allBadTimes{fileNum,iDate},estChInCortexA{iDate}(iRun,:),...
                estChInCortexB{iDate}(iRun,:),gammaRange,lowFreqRange);

            dataLen = size(amplitudeA,2);
            nChan      = min(size(amplitudeA,3),size(amplitudeB,3));

            amplitudeA = amplitudeA(:,:,1:nChan);
            amplitudeB = amplitudeB(:,:,1:nChan);
            phaseA     = phaseA(:,:,1:nChan);
            phaseB     = phaseB(:,:,1:nChan);


            winStart = 1:100e3:(dataLen-100e3);
            winEnd   = winStart+100e3-1;
            nWin     = numel(winStart);

            modIdxAllA2BShuffleT = NaN(5,nWin,nHigh,nLow,nChan); % size is shift length x # time windows x high freq x low freq x # channels
            modIdxAllB2AShuffleT = NaN(5,nWin,nHigh,nLow,nChan);
            modIdxAllA2AShuffleT = NaN(5,nWin,nHigh,nLow,nChan);
            modIdxAllB2BShuffleT = NaN(5,nWin,nHigh,nLow,nChan);

            % Create artificial surrogates for the data
            % Shuffle phase in windows
            clear phaseAShuff phaseBShuff
            phaseAShuff = ones([5,size(phaseA)],'single');
            phaseBShuff = ones([5,size(phaseB)],'single');
            winSize = 10e3;

            for iRep = 1:5
                rng('shuffle');
                comb1 = randperm(round(dataLen/10e3));
                rowIdx = 1;
                for iL = 1:length(comb1)
                    clear win1
                    win1 = ((comb1(iL)-1)*winSize+1 : (comb1(iL)-1)*winSize+winSize);
                    win1(win1>dataLen) = [];
                    numWin1 = length(win1);
                    phaseAShuff(iRep,:,rowIdx:rowIdx+numWin1-1, :) = phaseA(:,win1, :);
                    phaseBShuff(iRep,:,rowIdx:rowIdx+numWin1-1, :) = phaseB(:,win1, :);
                    rowIdx = rowIdx + numWin1;
                end
            end

            % Get the PAC

            % Calculate windows
            winStart = 1:100e3:(dataLen-100e3);
            winEnd   = winStart+100e3-1;
            nWin     = numel(winStart);

            highFreqAconst = parallel.pool.Constant(amplitudeA);
            highFreqBconst = parallel.pool.Constant(amplitudeB);
            lowFreqAconst  = parallel.pool.Constant(phaseAShuff);
            lowFreqBconst  = parallel.pool.Constant(phaseBShuff);

            comb = combvec(1:nHigh, 1:nLow, 1:nWin, 1:5)';
            nComb = size(comb,1);
            clear modShuffleA2B modShuffleB2A modShuffleA2A modShuffleB2B
            tic;
            parfor iC = 1:nComb
                iHigh   = comb(iC,1);
                iLow    = comb(iC,2);
                iWin    = comb(iC,3);
                iRep    = comb(iC,4);

                idx = winStart(iWin):winEnd(iWin);

                lowA  = lowFreqAconst.Value;
                lowB  = lowFreqBconst.Value;
                highA = highFreqAconst.Value;
                highB = highFreqBconst.Value;

                modShuffleA2B{iC} = getPhaseAmpCoupling(squeeze(lowB(iRep,iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)));
                modShuffleB2A{iC} = getPhaseAmpCoupling(squeeze(lowA(iRep,iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)));
                modShuffleA2A{iC} = getPhaseAmpCoupling(squeeze(lowA(iRep,iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)));
                modShuffleB2B{iC} = getPhaseAmpCoupling(squeeze(lowB(iRep,iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)));
            end
            toc;

            for iC = 1:nComb
                iHigh   = comb(iC,1);
                iLow    = comb(iC,2);
                iWin    = comb(iC,3);
                iRep    = comb(iC,4);

                modIdxAllA2BShuffleT(iRep,iWin,iHigh,iLow,:) = modShuffleA2B{iC};
                modIdxAllB2AShuffleT(iRep,iWin,iHigh,iLow,:) = modShuffleB2A{iC};
                modIdxAllA2AShuffleT(iRep,iWin,iHigh,iLow,:) = modShuffleA2A{iC};
                modIdxAllB2BShuffleT(iRep,iWin,iHigh,iLow,:) = modShuffleB2B{iC};
            end

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramCtrl_' num2str(fileNum) '.mat'],...
                'modIdxAllA2BShuffleT','modIdxAllB2AShuffleT','modIdxAllA2AShuffleT','modIdxAllB2BShuffleT');


            % tic;
            % highFreqAconst = parallel.pool.Constant(amplitudeA);
            % highFreqBconst = parallel.pool.Constant(amplitudeB);
            % for iW = 1: nWin
            %
            %     % winStart = 1:winSize(iW):(dataLen-winSize(iW));
            %     % winEnd   = winStart+winSize(iW)-1;
            %     % nWin     = min(10,numel(winStart));
            %
            %     % lowFreqAconst  = parallel.pool.Constant(phaseA);
            %     % lowFreqBconst  = parallel.pool.Constant(phaseB);
            %
            %     % Shuffle method 1: Segment the data into different 100 s windows,
            %     % shuffle the windows, and get the modulation index
            %     for iRep = 1:10
            %         disp(['Window length: ' num2str(winSize(iW)) ' Rep: ' num2str(iRep)]);
            %         rng('shuffle');
            %         comb1 = randperm(round(dataLen/winSize(iW)));
            %         clear newPhase gammaNew gammaNewFFT magGammaNew
            %
            %         % Shuffle phase in windows
            %         newAmpA = ones(size(amplitudeA));
            %         newAmpB = ones(size(amplitudeA));
            %         rowIdx = 1;
            %         for iL = 1:length(comb1)
            %             clear win1
            %             win1 = ((comb1(iL)-1)*winSize(iW)+1 : (comb1(iL)-1)*winSize(iW)+winSize(iW));
            %             win1(win1>dataLen) = [];
            %             numWin1 = length(win1);
            %             newAmpA(:,rowIdx:rowIdx+numWin1-1, :) = amplitudeA(:,win1, :);
            %             newAmpB(:,rowIdx:rowIdx+numWin1-1, :) = amplitudeB(:,win1, :);
            %             rowIdx = rowIdx + numWin1;
            %         end
            %
            %
            %         winStartConst  = parallel.pool.Constant(winStart);
            %         winEndConst    = parallel.pool.Constant(winEnd);
            %
            %         comb = combvec(1:nHigh, 1:nLow,1:nWin)';
            %         nComb = size(comb,1);
            %
            %         tic;
            %         parfor iShuffle = 1: nComb
            %             iHigh = comb(iShuffle,1);
            %             iLow  = comb(iShuffle,2);
            %             iWin  = comb(iShuffle,3);
            %
            %             lowA  = lowFreqAconst.Value;
            %             lowB  = lowFreqBconst.Value;
            %             highA = highFreqAconst.Value;
            %             highB = highFreqBconst.Value;
            %
            %             startVal = winStartConst.Value;
            %             endVal   = winEndConst.Value;
            %
            %             idx = startVal(iWin):endVal(iWin);
            %
            %             modShuffleA2B{iShuffle}  = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,:)),squeeze(highA(iHigh,idx,:)));
            %             modShuffleB2A{iShuffle}  = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,:)),squeeze(highB(iHigh,idx,:)));
            %             modShuffleA2A{iShuffle}  = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,:)),squeeze(highA(iHigh,idx,:)));
            %             modShuffleB2B{iShuffle}  = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,:)),squeeze(highB(iHigh,idx,:)));
            %
            %         end
            %         toc;
            %
            %         for iC = 1:nComb
            %             iHigh = comb(iC,1);
            %             iLow  = comb(iC,2);
            %             iWin  = comb(iC,3);
            %
            %             modIdxAllA2BShuffleT(iW,iRep,iWin,iHigh,iLow,:) = modShuffleA2B{iC};
            %             modIdxAllB2AShuffleT(iW,iRep,iWin,iHigh,iLow,:) = modShuffleB2A{iC};
            %             modIdxAllA2AShuffleT(iW,iRep,iWin,iHigh,iLow,:) = modShuffleA2A{iC};
            %             modIdxAllB2BShuffleT(iW,iRep,iWin,iHigh,iLow,:) = modShuffleB2B{iC};
            %         end
            %     end
            %
            % end

            % Circshift phase
            % highFreqAconst = parallel.pool.Constant(amplitudeA);
            % highFreqBconst = parallel.pool.Constant(amplitudeB);
            % for iWin = 1:nWin
            %     % Divide into epochs
            %     clear idx
            %     idx = winStart(iWin):winEnd(iWin);
            %     for iShift = 1:nShift
            %         clear newPhaseA newPhaseB
            %         % Shift the phase and calculate MI
            %         newPhaseA = circshift(phaseA,shiftLen(iShift),2);
            %         newPhaseB = circshift(phaseB,shiftLen(iShift),2);
            %
            %         lowFreqAconst  = parallel.pool.Constant(newPhaseA);
            %         lowFreqBconst  = parallel.pool.Constant(newPhaseB);
            %
            %         comb = combvec(1:nHigh, 1:nLow)';
            %         nComb = size(comb,1);
            %
            %         parfor iC = 1:nComb
            %             iHigh = comb(iC,1);
            %             iLow  = comb(iC,2);
            %
            %             lowA  = lowFreqAconst.Value;
            %             lowB  = lowFreqBconst.Value;
            %             highA = highFreqAconst.Value;
            %             highB = highFreqBconst.Value;
            %
            %             modShuffleA2B{iC}  = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,:)),squeeze(highA(iHigh,idx,:)));
            %             modShuffleB2A{iC}  = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,:)),squeeze(highB(iHigh,idx,:)));
            %             modShuffleA2A{iC}  = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,:)),squeeze(highA(iHigh,idx,:)));
            %             modShuffleB2B{iC}  = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,:)),squeeze(highB(iHigh,idx,:)));
            %         end
            %
            %         for iC = 1:nComb
            %             iHigh = comb(iC,1);
            %             iLow  = comb(iC,2);
            %
            %             modIdxAllA2BShuffleT(iShift,iWin,iHigh,iLow,:) = modShuffleA2B{iC};
            %             modIdxAllB2AShuffleT(iShift,iWin,iHigh,iLow,:) = modShuffleB2A{iC};
            %             modIdxAllA2AShuffleT(iShift,iWin,iHigh,iLow,:) = modShuffleA2A{iC};
            %             modIdxAllB2BShuffleT(iShift,iWin,iHigh,iLow,:) = modShuffleB2B{iC};
            %         end
            %
            %
            %     end
            %
            % end
            modIdxAllA2BShuffle{iRun,iDate} = modIdxAllA2BShuffleT;
            modIdxAllB2AShuffle{iRun,iDate} = modIdxAllB2AShuffleT;
            modIdxAllA2AShuffle{iRun,iDate} = modIdxAllA2AShuffleT;
            modIdxAllB2BShuffle{iRun,iDate} = modIdxAllB2BShuffleT;
            % toc;

        else
            clear allCtrlVals
            allCtrlVals = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramCtrl_' num2str(fileNum) '.mat']);
            modIdxAllA2BShuffle{iRun,iDate} = allCtrlVals.modIdxAllA2BShuffleT;
            modIdxAllB2AShuffle{iRun,iDate} = allCtrlVals.modIdxAllB2AShuffleT;
            modIdxAllA2AShuffle{iRun,iDate} = allCtrlVals.modIdxAllA2AShuffleT;
            modIdxAllB2BShuffle{iRun,iDate} = allCtrlVals.modIdxAllB2BShuffleT;
        end
    end
end
%%
modIdxAllA2B95 = reshape(cellfun(@(x) reshape(squeeze(median(x,5,'omitnan')),[size(x,1)*size(x,2) nHigh*nLow]),modIdxAllA2BShuffle,'UniformOutput',false),[1 numel(modIdxAllA2BShuffle)]);
zeroVals = cell2mat(cellfun(@(x) isempty(x), modIdxAllA2B95,'un',0));
modIdxAllA2B95(cell2mat(cellfun(@(x) isempty(x), modIdxAllA2B95,'un',0))) = [];
% modIdxAllA2B95 = cellfun(@(x) single(reshape(prctile(x,95,1),[nHigh nLow])),modIdxAllA2B95,'UniformOutput',0);

modIdxAllB2A95 = cellfun(@(x) reshape(squeeze(median(x,5,'omitnan')),[size(x,1)*size(x,2) nHigh*nLow]),modIdxAllB2AShuffle,'UniformOutput',false);
modIdxAllB2A95(cell2mat(cellfun(@(x) isempty(x), modIdxAllB2A95,'un',0))) = [];

% modIdxAllB2A95 = cellfun(@(x) single(reshape(prctile(x,95,1),[nHigh nLow])),modIdxAllB2A95,'UniformOutput',0);

modIdxAllA2A95 = cellfun(@(x) reshape(squeeze(median(x,5,'omitnan')),[size(x,1)*size(x,2) nHigh*nLow]),modIdxAllA2AShuffle,'UniformOutput',false);
modIdxAllA2A95(cell2mat(cellfun(@(x) isempty(x), modIdxAllA2A95,'un',0))) = [];

% modIdxAllA2A95 = cellfun(@(x) single(reshape(prctile(x,95,1),[nHigh nLow])),modIdxAllA2A95,'UniformOutput',0);

modIdxAllB2B95 = cellfun(@(x) reshape(squeeze(median(x,5,'omitnan')),[size(x,1)*size(x,2) nHigh*nLow]),modIdxAllB2BShuffle,'UniformOutput',false);
modIdxAllB2B95(cell2mat(cellfun(@(x) isempty(x), modIdxAllB2B95,'un',0))) = [];

% modIdxAllB2B95 = cellfun(@(x) single(reshape(prctile(x,95,1),[nHigh nLow])),modIdxAllB2B95,'UniformOutput',0);

% correctedA2B = cellfun(@(x,y) (squeeze(median(y,[1 4],'omitnan'))>x),modIdxAllA2B95,modIdxAllA2B,'un',0);
% correctedB2A = cellfun(@(x,y) (squeeze(median(y,[1 4],'omitnan'))>x),modIdxAllB2A95,modIdxAllB2A,'un',0);
% correctedA2A = cellfun(@(x,y) (squeeze(median(y,[1 4],'omitnan'))>x),modIdxAllA2A95,modIdxAllA2A,'un',0);
% correctedB2B = cellfun(@(x,y) (squeeze(median(y,[1 4],'omitnan'))>x),modIdxAllB2B95,modIdxAllB2B,'un',0);
modIdxAllA2BR = reshape(modIdxAllA2B,[1 numel(modIdxAllA2B)]); modIdxAllA2BR(cell2mat(cellfun(@(x) isempty(x), modIdxAllA2BR,'un',0))) = [];
modIdxAllB2AR = reshape(modIdxAllA2B,[1 numel(modIdxAllA2B)]); modIdxAllB2AR(cell2mat(cellfun(@(x) isempty(x), modIdxAllB2AR,'un',0))) = [];
modIdxAllA2AR = reshape(modIdxAllA2B,[1 numel(modIdxAllA2B)]); modIdxAllA2AR(cell2mat(cellfun(@(x) isempty(x), modIdxAllA2AR,'un',0))) = [];
modIdxAllB2BR = reshape(modIdxAllA2B,[1 numel(modIdxAllA2B)]); modIdxAllB2BR(cell2mat(cellfun(@(x) isempty(x), modIdxAllB2BR,'un',0))) = [];

%%
correctedA2B = cell2mat(cellfun(@(x,y) (reshape(squeeze(median(y,[1 4],'omitnan')),[1 nHigh*nLow])-mean(x,1,'omitnan'))./std(x),modIdxAllA2B95,modIdxAllA2BR,'un',0)');
correctedB2A = cell2mat(cellfun(@(x,y) (reshape(squeeze(median(y,[1 4],'omitnan')),[1 nHigh*nLow])-mean(x,1,'omitnan'))./std(x),modIdxAllB2A95,modIdxAllB2AR,'un',0)');
correctedA2A = cell2mat(cellfun(@(x,y) (reshape(squeeze(median(y,[1 4],'omitnan')),[1 nHigh*nLow])-mean(x,1,'omitnan'))./std(x),modIdxAllA2A95,modIdxAllA2AR,'un',0)');
correctedB2B = cell2mat(cellfun(@(x,y) (reshape(squeeze(median(y,[1 4],'omitnan')),[1 nHigh*nLow])-mean(x,1,'omitnan'))./std(x),modIdxAllB2B95,modIdxAllB2BR,'un',0)');

meanA2B = mean(correctedA2B,1); sigA2B = meanA2B>1.96;
meanB2A = mean(correctedB2A,1); sigB2A = meanB2A>1.96;
meanA2A = mean(correctedA2A,1); sigA2A = meanA2A>1.96;
meanB2B = mean(correctedB2B,1); sigB2B = meanB2B>1.96; 

commonPairs = sigA2B & sigB2A & sigA2A & sigB2A;

a2bAll = correctedA2B; a2bAll(:,~commonPairs) = NaN;
b2aAll = correctedB2A; b2aAll(:,~commonPairs) = NaN;
a2aAll = correctedA2A; a2aAll(:,~commonPairs) = NaN;
b2bAll = correctedB2B; b2bAll(:,~commonPairs) = NaN;

colIdx = [12:16 23:27];
figure;boxplot([median(a2bAll(:,colIdx),2,'omitnan') median(b2aAll(:,colIdx),2,'omitnan') ...
    median(a2aAll(:,colIdx),2,'omitnan') median(b2bAll(:,colIdx),2,'omitnan')],...
    {'A-B','B-A','A-A','B-B'});


figure; 
boxplot([median(correctedA2B(:,sigA2B),2) median(correctedB2A(:,sigB2A),2) ...
    median(correctedA2A(:,sigA2A),2) median(correctedB2B(:,sigB2B),2)],...
    {'A-B','B-A','A-A','B-B'});

connValsT = connValsAll(:,1:5);
connValsT = reshape(connValsT,[1 numel(connValsT)]);connValsT(zeroVals) = [];

distValsT = distSitesAll(:,1:5); 
distValsT = reshape(distValsT,[1 numel(distValsT)]); distValsT(zeroVals)= [];

[yIndex,edges] = discretize(connValsT,-0.4:0.1:0.8);
loc1 = yIndex<=8; loc2 = yIndex>8;
unlikePairs = (median(correctedA2B(:,sigA2B),2)+ median(correctedB2A(:,sigB2A),2))./2;
likePairs   = [median(correctedA2A(:,sigA2A),2); median(correctedB2B(:,sigB2B),2)];

figure;
subplot(121); histogram(unlikePairs(loc1),0:1:ceil(max(unlikePairs))); hold on;box off;
histogram(unlikePairs(loc2),0:1:ceil(max(unlikePairs))); hold on;box off;
% xlim([0 1e-4]); ylim([0 15]);
legend('FC<=0.3','FC>0.3','Location','northeast');

subplot(122); histogram(likePairs([loc1;loc1]),0:1:ceil(max(unlikePairs))); hold on;box off;
histogram(likePairs([loc2;loc2]),0:1:ceil(max(unlikePairs))); hold on;box off;
% xlim([0 1e-4]); ylim([0 15]);box off;
legend('FC<=0.3','FC>0.3','Location','northeast');

%%
figure;
subplot(221); scatter(connValsT,unlikePairs,'filled');
showLinearFit(connValsT',unlikePairs,-0.2,25,20); axis square;
xlabel('Functional connectivity'); ylabel('Z-scored Modulation index');
ylim([-2 35]); title('Between probes');

subplot(222); scatter([connValsT connValsT],likePairs,'filled'); 
showLinearFit([connValsT connValsT]',likePairs,-0.2,25,20); axis square;
xlabel('Functional connectivity'); ylabel('Z-scored Modulation index');
ylim([-2 35]); title('Within probe');

subplot(223); scatter(distValsT,unlikePairs,'filled');
showExpFit(distValsT',unlikePairs,15,25,20); axis square;
xlabel('Distance'); ylabel('Z-scored Modulation index');
ylim([-2 35]);

subplot(224); scatter([distValsT distValsT],likePairs,'filled');
showExpFit([distValsT distValsT]',likePairs,15,25,20); axis square;
xlabel('Distance'); ylabel('Z-scored Modulation index');
ylim([-2 35]);

[rhoBtw,pBtw] = partialcorr(connValsT',unlikePairs,distValsT');
[rhoLike,pLike] = partialcorr([connValsT connValsT]',likePairs,[distValsT distValsT]');

%% Flag the modulation indices which are higher than 5
highMIBtwFlag = unlikePairs>= 5; 

highA2BBtw = a2bAll(highMIBtwFlag,:);
highB2ABtw = b2aAll(highMIBtwFlag,:);

highA2ABtw = a2aAll(likePairs(1:55)>=5,:);
highB2BBtw = b2bAll(likePairs(56:end)>=5,:);

figure;
subplot(221); imagesc(lowFreqRange,gammaRange,squeeze(reshape(median(highA2BBtw,1,'omitnan'),[nHigh nLow])));
clim([2 8]); title('A-B'); ylabel('Amplitude'); xlabel('Phase'); axis square; set(gca,'YDir','normal');colorbar;

subplot(222); imagesc(lowFreqRange,gammaRange,squeeze(reshape(median(highB2ABtw,1,'omitnan'),[nHigh nLow])));
clim([2 8]); title('B-A'); ylabel('Amplitude'); xlabel('Phase'); axis square; set(gca,'YDir','normal');colorbar; 

subplot(223); imagesc(lowFreqRange,gammaRange,squeeze(reshape(median(highA2ABtw,1,'omitnan'),[nHigh nLow])));
clim([2 8]); title('A-A'); ylabel('Amplitude'); xlabel('Phase'); axis square; set(gca,'YDir','normal'); colorbar;

subplot(224); imagesc(lowFreqRange,gammaRange,squeeze(reshape(median(highB2BBtw,1,'omitnan'),[nHigh nLow])));
clim([2 8]); title('B-B'); ylabel('Amplitude'); xlabel('Phase'); axis square; set(gca,'YDir','normal'); colorbar;


%% 
modIdxAllA2BShuffleR = reshape(modIdxAllA2BShuffle,[1 numel(modIdxAllA2BShuffle)]); modIdxAllA2BShuffleR(cell2mat(cellfun(@(x) isempty(x), modIdxAllA2BShuffleR,'un',0))) = [];
modIdxAllA2BShuffleR = cellfun(@(x) reshape(x,[size(x,1)*size(x,2) nHigh*nLow size(x,5)]),modIdxAllA2BShuffleR,'un',0);

modIdxAllB2AShuffleR = reshape(modIdxAllB2AShuffle,[1 numel(modIdxAllB2AShuffle)]); modIdxAllB2AShuffleR(cell2mat(cellfun(@(x) isempty(x), modIdxAllB2AShuffleR,'un',0))) = [];
modIdxAllB2AShuffleR = cellfun(@(x) reshape(x,[size(x,1)*size(x,2) nHigh*nLow size(x,5)]),modIdxAllB2AShuffleR,'un',0);

modIdxAllA2AShuffleR = reshape(modIdxAllA2AShuffle,[1 numel(modIdxAllA2AShuffle)]); modIdxAllA2AShuffleR(cell2mat(cellfun(@(x) isempty(x), modIdxAllA2AShuffleR,'un',0))) = [];
modIdxAllA2AShuffleR = cellfun(@(x) reshape(x,[size(x,1)*size(x,2) nHigh*nLow size(x,5)]),modIdxAllA2AShuffleR,'un',0);

modIdxAllB2BShuffleR = reshape(modIdxAllB2BShuffle,[1 numel(modIdxAllB2BShuffle)]); modIdxAllB2BShuffleR(cell2mat(cellfun(@(x) isempty(x), modIdxAllB2BShuffleR,'un',0))) = [];
modIdxAllB2BShuffleR = cellfun(@(x) reshape(x,[size(x,1)*size(x,2) nHigh*nLow size(x,5)]),modIdxAllB2BShuffleR,'un',0);

a2bLaminar = (cellfun(@(x,y)(reshape(squeeze(mean(y,1,'omitnan')),[nHigh*nLow size(y,4) ])-(squeeze(mean(x,1,'omitnan'))))./(squeeze(std(x))),modIdxAllA2BShuffleR,modIdxAllA2BR,'un',0));
b2aLaminar = (cellfun(@(x,y)(reshape(squeeze(mean(y,1,'omitnan')),[nHigh*nLow size(y,4) ])-(squeeze(mean(x,1,'omitnan'))))./(squeeze(std(x))),modIdxAllB2AShuffleR,modIdxAllB2AR,'un',0));
a2aLaminar = (cellfun(@(x,y)(reshape(squeeze(mean(y,1,'omitnan')),[nHigh*nLow size(y,4) ])-(squeeze(mean(x,1,'omitnan'))))./(squeeze(std(x))),modIdxAllA2AShuffleR,modIdxAllA2AR,'un',0));
b2bLaminar = (cellfun(@(x,y)(reshape(squeeze(mean(y,1,'omitnan')),[nHigh*nLow size(y,4) ])-(squeeze(mean(x,1,'omitnan'))))./(squeeze(std(x))),modIdxAllB2BShuffleR,modIdxAllB2BR,'un',0));

%%
colVals = zeros(1,143); colVals(colIdx) =1; 
a2bLaminarSuper = cell2mat(cellfun(@(x) median(x(:,1:6),2,'omitnan'),a2bLaminar,'UniformOutput',0));
a2bLaminarSuper(~colVals,:) = NaN;   a2bLaminarSuper = reshape(a2bLaminarSuper,[nHigh nLow numel(a2bLaminar)]); 
a2bLaminarMid = cell2mat(cellfun(@(x) median(x(:,7:12),2,'omitnan'),a2bLaminar,'UniformOutput',0));
a2bLaminarMid (~colVals,:) = NaN;   a2bLaminarMid  = reshape(a2bLaminarMid ,[nHigh nLow numel(a2bLaminar)]); 
a2bLaminarDeep = cell2mat(cellfun(@(x) median(x(:,13:end),2,'omitnan'),a2bLaminar,'UniformOutput',0));
a2bLaminarDeep(~colVals,:) = NaN;   a2bLaminarDeep = reshape(a2bLaminarDeep,[nHigh nLow numel(a2bLaminar)]); 

b2aLaminarSuper = cell2mat(cellfun(@(x) median(x(:,1:6),2,'omitnan'),b2aLaminar,'UniformOutput',0));
b2aLaminarSuper(~colVals,:) = NaN;   b2aLaminarSuper = reshape(b2aLaminarSuper,[nHigh nLow numel(a2bLaminar)]); 
b2aLaminarMid = cell2mat(cellfun(@(x) median(x(:,7:12),2,'omitnan'),b2aLaminar,'UniformOutput',0));
b2aLaminarMid (~colVals,:) = NaN;   b2aLaminarMid  = reshape(b2aLaminarMid ,[nHigh nLow numel(a2bLaminar)]); 
b2aLaminarDeep = cell2mat(cellfun(@(x) median(x(:,13:end),2,'omitnan'),b2aLaminar,'UniformOutput',0));
b2aLaminarDeep(~colVals,:) = NaN;   b2aLaminarDeep = reshape(b2aLaminarDeep,[nHigh nLow numel(a2bLaminar)]); 

a2aLaminarSuper = cell2mat(cellfun(@(x) median(x(:,1:6),2,'omitnan'),a2aLaminar,'UniformOutput',0));
a2aLaminarSuper(~colVals,:) = NaN;   a2aLaminarSuper = reshape(a2aLaminarSuper,[nHigh nLow numel(a2bLaminar)]); 
a2aLaminarMid = cell2mat(cellfun(@(x) median(x(:,7:12),2,'omitnan'),a2aLaminar,'UniformOutput',0));
a2aLaminarMid (~colVals,:) = NaN;   a2aLaminarMid  = reshape(a2aLaminarMid ,[nHigh nLow numel(a2bLaminar)]); 
a2aLaminarDeep = cell2mat(cellfun(@(x) median(x(:,13:end),2,'omitnan'),a2aLaminar,'UniformOutput',0));
a2aLaminarDeep(~colVals,:) = NaN;   a2aLaminarDeep = reshape(a2aLaminarDeep,[nHigh nLow numel(a2bLaminar)]); 

b2bLaminarSuper = cell2mat(cellfun(@(x) median(x(:,1:6),2,'omitnan'),b2bLaminar,'UniformOutput',0));
b2bLaminarSuper(~colVals,:) = NaN;   b2bLaminarSuper = reshape(b2bLaminarSuper,[nHigh nLow numel(a2bLaminar)]); 
b2bLaminarMid = cell2mat(cellfun(@(x) median(x(:,7:12),2,'omitnan'),b2bLaminar,'UniformOutput',0));
b2bLaminarMid(~colVals,:) = NaN;   b2bLaminarMid = reshape(b2bLaminarMid,[nHigh nLow numel(a2bLaminar)]); 
b2bLaminarDeep = cell2mat(cellfun(@(x) median(x(:,13:end),2,'omitnan'),b2bLaminar,'UniformOutput',0));
b2bLaminarDeep(~colVals,:) = NaN;   b2bLaminarDeep = reshape(b2bLaminarDeep,[nHigh nLow numel(a2bLaminar)]); 

%% 
figure; 
subplot(341); imagesc(lowFreqRange,gammaRange,squeeze(median(a2bLaminarSuper,3,'omitnan'))); clim([0 5]); title('A-B'); axis square; set(gca,'YDir','normal'); ylabel('Superficial');
subplot(342); imagesc(lowFreqRange,gammaRange,squeeze(median(b2aLaminarSuper,3,'omitnan'))); clim([0 5]); title('B-A');axis square; set(gca,'YDir','normal'); 
subplot(343); imagesc(lowFreqRange,gammaRange,squeeze(median(a2aLaminarSuper,3,'omitnan'))); clim([0 5]); title('A-A');axis square; set(gca,'YDir','normal'); 
subplot(344); imagesc(lowFreqRange,gammaRange,squeeze(median(b2bLaminarSuper,3,'omitnan'))); clim([0 5]); title('B-B');colorbar;axis square; set(gca,'YDir','normal'); 

subplot(345); imagesc(lowFreqRange,gammaRange,squeeze(median(a2bLaminarMid,3,'omitnan'))); clim([0 5]); title('A-B'); ylabel('Middle');axis square; set(gca,'YDir','normal'); 
subplot(346); imagesc(lowFreqRange,gammaRange,squeeze(median(b2aLaminarMid,3,'omitnan'))); clim([0 5]); title('B-A');axis square; set(gca,'YDir','normal'); 
subplot(347); imagesc(lowFreqRange,gammaRange,squeeze(median(a2aLaminarMid,3,'omitnan'))); clim([0 5]); title('A-A');axis square; set(gca,'YDir','normal'); 
subplot(348); imagesc(lowFreqRange,gammaRange,squeeze(median(b2bLaminarMid,3,'omitnan'))); clim([0 5]); title('B-B');colorbar;axis square; set(gca,'YDir','normal'); 

subplot(349); imagesc(lowFreqRange,gammaRange,squeeze(median(a2bLaminarDeep,3,'omitnan'))); clim([0 5]); title('A-B'); ylabel('Deep');axis square; set(gca,'YDir','normal'); 
subplot(3,4,10); imagesc(lowFreqRange,gammaRange,squeeze(median(b2aLaminarDeep,3,'omitnan'))); clim([0 5]); title('B-A');axis square; set(gca,'YDir','normal'); 
subplot(3,4,11); imagesc(lowFreqRange,gammaRange,squeeze(median(a2aLaminarDeep,3,'omitnan'))); clim([0 5]); title('A-A');axis square; set(gca,'YDir','normal'); 
subplot(3,4,12); imagesc(lowFreqRange,gammaRange,squeeze(median(b2bLaminarDeep,3,'omitnan'))); clim([0 5]); title('B-B');colorbar;axis square; set(gca,'YDir','normal'); 
xlabel('Phase'); ylabel('Amplitude');


%% Plot
figure; subplot(131); imagesc(lowFreqRange,gammaRange,reshape(median(a2bLaminar{1}(:,1:6),2,'omitnan'),[11 13]));
clim([1.96 5]); colorbar;
subplot(132); imagesc(lowFreqRange,gammaRange,reshape(median(a2bLaminar{1}(:,7:12),2,'omitnan'),[11 13]));
clim([1.96 5]);colorbar;
subplot(133); imagesc(lowFreqRange,gammaRange,reshape(median(a2bLaminar{1}(:,13:end),2,'omitnan'),[11 13]));
clim([1.96 5]);colorbar;

%%
a2bSuper = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,1:6),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,1:6),[1 3],'omitnan')))./squeeze(std(x(:,:,1:6),[],[1 3])),modIdxAllA2BShuffleR,modIdxAllA2BR,'un',0)'); 
a2bSuper(:,~colVals) = NaN;  a2bSuper = reshape(a2bSuper,[numel(a2bLaminar) nHigh nLow]); 
a2bMid   = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,7:12),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,7:12),[1 3],'omitnan')))./squeeze(std(x(:,:,7:12),[],[1 3])),modIdxAllA2BShuffleR,modIdxAllA2BR,'un',0)');
a2bMid(:,~colVals) = NaN;  a2bMid = reshape(a2bMid,[numel(a2bLaminar) nHigh nLow]); 
a2bDeep  = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,13:end),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,13:end),[1 3],'omitnan')))./squeeze(std(x(:,:,13:end),[],[1 3])),modIdxAllA2BShuffleR,modIdxAllA2BR,'un',0)');
a2bDeep(:,~colVals) = NaN;  a2bDeep = reshape(a2bDeep,[numel(a2bLaminar) nHigh nLow]); 

b2aSuper = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,1:6),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,1:6),[1 3],'omitnan')))./squeeze(std(x(:,:,1:6),[],[1 3])),modIdxAllB2AShuffleR,modIdxAllB2AR,'un',0)');
b2aSuper(:,~colVals) = NaN;  b2aSuper = reshape(b2aSuper,[numel(a2bLaminar) nHigh nLow]); 
b2aMid   = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,7:12),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,7:12),[1 3],'omitnan')))./squeeze(std(x(:,:,7:12),[],[1 3])),modIdxAllB2AShuffleR,modIdxAllB2AR,'un',0)');
b2aMid (:,~colVals) = NaN;  b2aMid  = reshape(b2aMid ,[numel(a2bLaminar) nHigh nLow]); 
b2aDeep  = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,13:end),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,13:end),[1 3],'omitnan')))./squeeze(std(x(:,:,13:end),[],[1 3])),modIdxAllB2AShuffleR,modIdxAllB2AR,'un',0)');
b2aDeep(:,~colVals) = NaN;  b2aDeep = reshape(b2aDeep,[numel(a2bLaminar) nHigh nLow]); 

a2aSuper = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,1:6),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,1:6),[1 3],'omitnan')))./squeeze(std(x(:,:,1:6),[],[1 3])),modIdxAllA2AShuffleR,modIdxAllA2AR,'un',0)');
a2aSuper(:,~colVals) = NaN;  a2aSuper = reshape(a2aSuper,[numel(a2bLaminar) nHigh nLow]); 
a2aMid   = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,7:12),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,7:12),[1 3],'omitnan')))./squeeze(std(x(:,:,7:12),[],[1 3])),modIdxAllA2AShuffleR,modIdxAllA2AR,'un',0)');
a2aMid (:,~colVals) = NaN;  a2aMid  = reshape(a2aMid ,[numel(a2bLaminar) nHigh nLow]); 
a2aDeep  = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,13:end),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,13:end),[1 3],'omitnan')))./squeeze(std(x(:,:,13:end),[],[1 3])),modIdxAllA2AShuffleR,modIdxAllA2AR,'un',0)');
a2aDeep(:,~colVals) = NaN;  a2aDeep = reshape(a2aDeep,[numel(a2bLaminar) nHigh nLow]); 

b2bSuper = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,1:6),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,1:6),[1 3],'omitnan')))./squeeze(std(x(:,:,1:6),[],[1 3])),modIdxAllB2BShuffleR,modIdxAllB2BR,'un',0)');
b2bSuper(:,~colVals) = NaN;  b2bSuper = reshape(b2bSuper,[numel(a2bLaminar) nHigh nLow]); 
b2bMid   = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,7:12),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,7:12),[1 3],'omitnan')))./squeeze(std(x(:,:,7:12),[],[1 3])),modIdxAllB2BShuffleR,modIdxAllB2BR,'un',0)');
b2bMid (:,~colVals) = NaN;  b2bMid = reshape(b2bMid ,[numel(a2bLaminar) nHigh nLow]); 
b2bDeep  = cell2mat(cellfun(@(x,y)(reshape(squeeze(mean(y(:,:,:,13:end),[1 4],'omitnan')),[nHigh*nLow,1])' - squeeze(mean(x(:,:,13:end),[1 3],'omitnan')))./squeeze(std(x(:,:,13:end),[],[1 3])),modIdxAllB2BShuffleR,modIdxAllB2BR,'un',0)');
b2bDeep(:,~colVals) = NaN;  b2bDeep = reshape(b2bDeep,[numel(a2bLaminar) nHigh nLow]); 

%% 
figure; 
subplot(341); imagesc(lowFreqRange,gammaRange,squeeze(median(a2bSuper,1,'omitnan'))); clim([0 5]); title('A-B'); axis square; set(gca,'YDir','normal'); ylabel('Superficial');
subplot(342); imagesc(lowFreqRange,gammaRange,squeeze(median(b2aSuper,1,'omitnan'))); clim([0 5]); title('B-A');axis square; set(gca,'YDir','normal'); 
subplot(343); imagesc(lowFreqRange,gammaRange,squeeze(median(a2aSuper,1,'omitnan'))); clim([0 5]); title('A-A');axis square; set(gca,'YDir','normal'); 
subplot(344); imagesc(lowFreqRange,gammaRange,squeeze(median(b2bSuper,1,'omitnan'))); clim([0 5]); title('B-B');colorbar;axis square; set(gca,'YDir','normal'); 

subplot(345); imagesc(lowFreqRange,gammaRange,squeeze(median(a2bMid,1,'omitnan'))); clim([0 5]); title('A-B'); ylabel('Middle');axis square; set(gca,'YDir','normal'); 
subplot(346); imagesc(lowFreqRange,gammaRange,squeeze(median(b2aMid,1,'omitnan'))); clim([0 5]); title('B-A');axis square; set(gca,'YDir','normal'); 
subplot(347); imagesc(lowFreqRange,gammaRange,squeeze(median(a2aMid,1,'omitnan'))); clim([0 5]); title('A-A');axis square; set(gca,'YDir','normal'); 
subplot(348); imagesc(lowFreqRange,gammaRange,squeeze(median(b2bMid,1,'omitnan'))); clim([0 5]); title('B-B');colorbar;axis square; set(gca,'YDir','normal'); 

subplot(349); imagesc(lowFreqRange,gammaRange,squeeze(median(a2bDeep,1,'omitnan'))); clim([0 5]); title('A-B'); ylabel('Deep');axis square; set(gca,'YDir','normal'); 
subplot(3,4,10); imagesc(lowFreqRange,gammaRange,squeeze(median(b2aDeep,1,'omitnan'))); clim([0 5]); title('B-A');axis square; set(gca,'YDir','normal'); 
subplot(3,4,11); imagesc(lowFreqRange,gammaRange,squeeze(median(a2aDeep,1,'omitnan'))); clim([0 5]); title('A-A');axis square; set(gca,'YDir','normal'); 
subplot(3,4,12); imagesc(lowFreqRange,gammaRange,squeeze(median(b2bDeep,1,'omitnan'))); clim([0 5]); title('B-B');colorbar;axis square; set(gca,'YDir','normal'); 
xlabel('Phase'); ylabel('Amplitude');


%%
figure;
for iPlot = 1:3
    switch iPlot
        case 1
            a2b = cell2mat(cellfun(@ (x) median(x(sigA2B,1:6),2,'omitnan'), a2bLaminar,'un',0));
            b2a = cell2mat(cellfun(@ (x) median(x(sigB2A,1:6),2,'omitnan'), b2aLaminar,'un',0));
            a2a = cell2mat(cellfun(@ (x) median(x(sigA2A,1:6),2,'omitnan'), a2aLaminar,'un',0));
            b2b = cell2mat(cellfun(@ (x) median(x(sigB2B,1:6),2,'omitnan'), b2bLaminar,'un',0));
            pltTitle = 'Superficial';
        case 2
            a2b = cell2mat(cellfun(@ (x) median(x(sigA2B,7:12), 2,'omitnan'), a2bLaminar,'un',0));
            b2a = cell2mat(cellfun(@ (x) median(x(sigB2A,7:12), 2,'omitnan'), b2aLaminar,'un',0));
            a2a = cell2mat(cellfun(@ (x) median(x(sigA2A,7:12), 2,'omitnan'), a2aLaminar,'un',0));
            b2b = cell2mat(cellfun(@ (x) median(x(sigB2B,7:12), 2,'omitnan'), b2bLaminar,'un',0));
            pltTitle = 'Middle';
        case 3
            a2b = cell2mat(cellfun(@ (x) median(x(sigA2B,13:end), 2,'omitnan'), a2bLaminar,'un',0));
            b2a = cell2mat(cellfun(@ (x) median(x(sigB2A,13:end), 2,'omitnan'), b2aLaminar,'un',0));
            a2a = cell2mat(cellfun(@ (x) median(x(sigA2A,13:end), 2,'omitnan'), a2aLaminar,'un',0));
            b2b = cell2mat(cellfun(@ (x) median(x(sigB2B,13:end), 2,'omitnan'), b2bLaminar,'un',0));
            pltTitle = 'Deep';
    end


    % figure;
    % subplot(221); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
    % subplot(222); showLinearFit(connValsT,squeeze(b2a(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
    % subplot(223); showLinearFit(connValsT,squeeze(a2a(iLayer,:)),0.3,1.85e-4,1.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('A--A'); axis square;
    % subplot(224); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('B--B'); axis square;
    % xlabel('Functional connectivity'); ylabel('Modulation index')
    % sgtitle([pltTitle '- ' layerTitle]);

    subplot(1,3,iPlot)
    imagesc([[mean(a2b,2);NaN(3,1)]  [mean(b2a,2);NaN]  [mean(a2a,2)]  [mean(b2b,2);NaN(3,1)]]);%,{'A-B','B-A','A-A','B-B'});
    title(pltTitle ); box off; clim([1.96 5]);%ylim([0 2e-4]);
    xticks(1:4); xticklabels({'A-B','B-A','A-A','B-B'});
end

%%
figure;
for iPlot = 1:3
    switch iPlot
        case 1
            a2b = cell2mat(cellfun(@ (x) median(x(sigA2B,1:6),[1 2],'omitnan'), a2bLaminar,'un',0)');
            b2a = cell2mat(cellfun(@ (x) median(x(sigB2A,1:6),[1 2],'omitnan'), b2aLaminar,'un',0)');
            a2a = cell2mat(cellfun(@ (x) median(x(sigA2A,1:6),[1 2],'omitnan'), a2aLaminar,'un',0)');
            b2b = cell2mat(cellfun(@ (x) median(x(sigB2B,1:6),[1 2],'omitnan'), b2bLaminar,'un',0)');
            pltTitle = 'Superficial';
        case 2
            a2b = cell2mat(cellfun(@ (x) median(x(sigA2B,7:12),[1 2],'omitnan'), a2bLaminar,'un',0)');
            b2a = cell2mat(cellfun(@ (x) median(x(sigB2A,7:12),[1 2],'omitnan'), b2aLaminar,'un',0)');
            a2a = cell2mat(cellfun(@ (x) median(x(sigA2A,7:12),[1 2],'omitnan'), a2aLaminar,'un',0)');
            b2b = cell2mat(cellfun(@ (x) median(x(sigB2B,7:12),[1 2],'omitnan'), b2bLaminar,'un',0)');
            pltTitle = 'Middle';
        case 3
            a2b = cell2mat(cellfun(@ (x) median(x(sigA2B,13:end),[1 2],'omitnan'), a2bLaminar,'un',0)');
            b2a = cell2mat(cellfun(@ (x) median(x(sigB2A,13:end),[1 2],'omitnan'), b2aLaminar,'un',0)');
            a2a = cell2mat(cellfun(@ (x) median(x(sigA2A,13:end),[1 2],'omitnan'), a2aLaminar,'un',0)');
            b2b = cell2mat(cellfun(@ (x) median(x(sigB2B,13:end),[1 2],'omitnan'), b2bLaminar,'un',0)');
            pltTitle = 'Deep';

    end


    % figure;
    % subplot(221); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
    % subplot(222); showLinearFit(connValsT,squeeze(b2a(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
    % subplot(223); showLinearFit(connValsT,squeeze(a2a(iLayer,:)),0.3,1.85e-4,1.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('A--A'); axis square;
    % subplot(224); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('B--B'); axis square;
    % xlabel('Functional connectivity'); ylabel('Modulation index')
    % sgtitle([pltTitle '- ' layerTitle]);

    subplot(1,3,iPlot)
    boxplot([a2b  b2a  a2a  b2b],{'A-B','B-A','A-A','B-B'});
    title(pltTitle ); box off; ylim([-1 4]);%ylim([0 2e-4]);
end


%% Getting the shuffled comodulogram - circshift method
winSize   = 100e3;
stepSize  = 100e3;
nHigh     = size(gammaRange,2);
nLow      = size(lowFreqRange,2);

for iDate = 2:3%1:size(allDates,1)
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        clear fileNum chA chB amplitudeA amplitudeB phaseA phaseB...
            highFreqAconst highFreqBconst lowFreqAconst lowFreqBconst...
            modCircleA2B modCircleB2A modCircleA2A modCircleB2B
        fileNum = datFileNum(iRun);

        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);
        if chA(1)== 0 || chB(1)==0; continue; end
        fileMat = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'];

        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file') %(~ismember('modIdxAllA2BCircleT', who('-file', fileMat)))

            clc; disp(['Getting control data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

            [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
                allProbeData{fileNum,iDate}.probe1Ch,allProbeData{fileNum,iDate}.probe2Ch,badElecA{fileNum,iDate},...
                badElecB{fileNum,iDate},allBadTimes{fileNum,iDate},estChInCortexA{iDate}(iRun,:),...
                estChInCortexB{iDate}(iRun,:),gammaRange,lowFreqRange);

            dataLen = size(amplitudeA,2);
            nChan      = min(size(amplitudeA,3),size(amplitudeB,3));
            amplitudeA = amplitudeA(:,:,1:nChan);
            amplitudeB = amplitudeB(:,:,1:nChan);
            phaseA     = phaseA(:,:,1:nChan);
            phaseB     = phaseB(:,:,1:nChan);

            % Calculate windows
            winStart = 1:stepSize:(dataLen-winSize);
            winEnd   = winStart+winSize-1;
            nWin     = numel(winStart);

            % % Shuffle method 1: Segment the data into different 100 s windows,
            % % shuffle the windows, and get the modulation index
            % winSize   = 100e3;
            % stepSize  = 100e3;
            %
            % % Generate all shuffle pairs (iWin, jWin) where iWin ~= jWin
            % [w1, w2]     = ndgrid(1:nWin, 1:nWin);
            % shufflePairs = [w1(:), w2(:)];
            % shufflePairs = shufflePairs(shufflePairs(:,1) ~= shufflePairs(:,2), :);
            % nShuffle     = size(shufflePairs,1);
            %
            % highFreqAconst = parallel.pool.Constant(amplitudeA);
            % highFreqBconst = parallel.pool.Constant(amplitudeB);
            % lowFreqAconst  = parallel.pool.Constant(phaseA);
            % lowFreqBconst  = parallel.pool.Constant(phaseB);
            %
            % comb = combvec(1:nHigh, 1:nLow, 1:nShuffle)';
            % nComb = size(comb,1);
            %
            % tic
            % parfor iShuffle = 1: nComb
            %     iHigh = comb(iShuffle,1);
            %     iLow  = comb(iShuffle,2);
            %     iRow  = comb(iShuffle,3);
            %
            %     idx1 = winStart(shufflePairs(iRow,1)):winEnd(shufflePairs(iRow,1));
            %     idx2 = winStart(shufflePairs(iRow,2)):winEnd(shufflePairs(iRow,2));
            %
            %     lowA  = lowFreqAconst.Value;
            %     lowB  = lowFreqBconst.Value;
            %     highA = highFreqAconst.Value;
            %     highB = highFreqBconst.Value;
            %
            %     modShuffleA2B{iShuffle} = getPhaseAmpCoupling(squeeze(lowB(iLow,idx1,:)),squeeze(highA(iHigh,idx2,:)));
            %     modShuffleB2A{iShuffle} = getPhaseAmpCoupling(squeeze(lowA(iLow,idx1,:)),squeeze(highB(iHigh,idx2,:)));
            %     modShuffleA2A{iShuffle} = getPhaseAmpCoupling(squeeze(lowA(iLow,idx1,:)),squeeze(highA(iHigh,idx2,:)));
            %     modShuffleB2B{iShuffle} = getPhaseAmpCoupling(squeeze(lowB(iLow,idx1,:)),squeeze(highB(iHigh,idx2,:)));
            % end
            % toc;
            %
            % for iC = 1:nComb
            %     iHigh = comb(iC,1);
            %     iLow  = comb(iC,2);
            %     iRow  = comb(iC,3);
            %
            %     modIdxAllA2BShuffleT(iRow,iHigh,iLow,:) = modShuffleA2B{iC};
            %     modIdxAllB2AShuffleT(iRow,iHigh,iLow,:) = modShuffleB2A{iC};
            %     modIdxAllA2AShuffleT(iRow,iHigh,iLow,:) = modShuffleA2A{iC};
            %     modIdxAllB2BShuffleT(iRow,iHigh,iLow,:) = modShuffleB2B{iC};
            % end
            %
            % modIdxAllA2BShuffle{iRun,iDate} = modIdxAllA2BShuffleT;
            % modIdxAllB2AShuffle{iRun,iDate} = modIdxAllB2AShuffleT;
            % modIdxAllA2AShuffle{iRun,iDate} = modIdxAllA2AShuffleT;
            % modIdxAllB2BShuffle{iRun,iDate} = modIdxAllB2BShuffleT;

            % Shuffle method 2: Circular shifting the data
            shiftLen = [0.01 0.1 0.2 0.3 0.4 0.5].*dataLen;
            nShift   = length(shiftLen);

            % Calculate the modulation index for a shifted distribution
            ampShiftA = NaN(nShift,nHigh,dataLen,nChan,'single');
            ampShiftB = NaN(nShift,nHigh,dataLen,nChan,'single');


            for iShift = 1:nShift
                ampShiftA(iShift,:,:,:) = single(circshift(amplitudeA,[0,shiftLen(iShift),0]));
                ampShiftB(iShift,:,:,:) = single(circshift(amplitudeB,[0,shiftLen(iShift),0]));
            end

            highFreqAconst = parallel.pool.Constant(ampShiftA);
            highFreqBconst = parallel.pool.Constant(ampShiftB);
            lowFreqAconst  = parallel.pool.Constant(phaseA);
            lowFreqBconst  = parallel.pool.Constant(phaseB);

            combShift  = combvec(1:nHigh, 1:nLow, 1:nShift,1:nWin)';
            nCombShift = size(combShift,1);

            tic;
            parfor iShuffle = 1:nCombShift
                iHigh   = combShift(iShuffle,1);
                iLow    = combShift(iShuffle,2);
                iShift  = combShift(iShuffle,3);
                iWin    = combShift(iShuffle,4);

                idx = winStart(iWin):winEnd(iWin);

                lowA  = lowFreqAconst.Value;
                lowB  = lowFreqBconst.Value;
                highA = highFreqAconst.Value;
                highB = highFreqBconst.Value;

                modCircleA2B{iShuffle} = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,:)),squeeze(highA(iShift,iHigh,idx,:)));
                modCircleB2A{iShuffle} = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,:)),squeeze(highB(iShift,iHigh,idx,:)));
                modCircleA2A{iShuffle} = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,:)),squeeze(highA(iShift,iHigh,idx,:)));
                modCircleB2B{iShuffle} = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,:)),squeeze(highB(iShift,iHigh,idx,:)));
            end
            toc;
            clear modIdxAllA2BCircleT modIdxAllB2ACircleT modIdxAllA2ACircleT modIdxAllB2BCircleT
            for iC = 1:nCombShift
                iHigh   = combShift(iC,1);
                iLow    = combShift(iC,2);
                iShift  = combShift(iC,3);
                iWin    = combShift(iC,4);

                modIdxAllA2BCircleT(iShift,iWin,iHigh,iLow,:) = modCircleA2B{iC};
                modIdxAllB2ACircleT(iShift,iWin,iHigh,iLow,:) = modCircleB2A{iC};
                modIdxAllA2ACircleT(iShift,iWin,iHigh,iLow,:) = modCircleA2A{iC};
                modIdxAllB2BCircleT(iShift,iWin,iHigh,iLow,:) = modCircleB2B{iC};
            end

            modIdxAllA2BCircle{iRun,iDate}  = modIdxAllA2BCircleT;
            modIdxAllB2ACircle{iRun,iDate} = modIdxAllB2ACircleT;
            modIdxAllA2ACircle{iRun,iDate} = modIdxAllA2ACircleT;
            modIdxAllB2BCircle{iRun,iDate} = modIdxAllB2BCircleT;

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
                'modIdxAllA2BCircleT','modIdxAllB2ACircleT','modIdxAllA2ACircleT','modIdxAllB2BCircleT','-append');
        else
            clear vars
            vars = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
            modIdxAllA2BCircle{iRun,iDate} = vars.modIdxAllA2BCircleT;
            modIdxAllB2ACircle{iRun,iDate} = vars.modIdxAllB2ACircleT;
            modIdxAllA2ACircle{iRun,iDate} = vars.modIdxAllA2ACircleT;
            modIdxAllB2BCircle{iRun,iDate} = vars.modIdxAllB2BCircleT;
        end
    end
    % 'modIdxAllA2BShuffleT','modIdxAllB2AShuffleT','modIdxAllA2AShuffleT','modIdxAllB2BShuffleT',...
end



%% Plot the comodulogram
avgModA2B = squeeze(mean(modIdxAllA2BShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
avgModB2A = squeeze(mean(modIdxAllB2AShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
figure; contourf(lowFreqRange,gammaRange,median(avgModA2B,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;
figure; contourf(lowFreqRange,gammaRange,median(avgModB2A,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;

for iPlot = 1:2
    switch iPlot
        case 1
            plotVar = avgModA2B;
            figTitle = 'A--B';
        case 2
            plotVar = avgModB2A;
            figTitle = 'B--A';
    end

    figure;
    for iCh = 1: size(plotVar,3)
        subplot(size(plotVar,3),1,iCh)
        contourf(lowFreqRange,gammaRange,plotVar(:,:,iCh),'lines','none'); colormap jet; clim([0 1e-4])
        if iCh~=size(plotVar,3)
            xticklabels({}); yticklabels({});
        end
    end
    sgtitle(figTitle);
end

% Plot for the control/shuffled distribution
% avgShuffledModA2B = squeeze(mean(modIdxAllA2BShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
% avgShuffledModB2A = squeeze(mean(modIdxAllB2AShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
% figure; contourf(lowFreqRange,gammaRange,median(avgShuffledModA2B,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;
% figure; contourf(lowFreqRange,gammaRange,median(avgShuffledModB2A,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;


% Plot for within probe
avgModA2A = squeeze(mean(modIdxAllA2AShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
avgModB2B = squeeze(mean(modIdxAllB2BShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
figure; contourf(lowFreqRange,gammaRange,median(avgModA2A,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;
figure; contourf(lowFreqRange,gammaRange,median(avgModB2B,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;

%% Plot the comodulogram -circularshift

for iPlot = 1:4
    switch iPlot
        case 1
            plotVar = modIdxAllA2BCircle{iRun,iDate};
            figTitle = 'A--B';
        case 2
            plotVar = modIdxAllB2ACircle{iRun,iDate};
            figTitle = 'B--A';
        case 3
            plotVar = modIdxAllA2ACircle{iRun,iDate};
            figTitle = 'A--A';
        case 4
            plotVar = modIdxAllB2BCircle{iRun,iDate};
            figTitle = 'B--B';
    end

    figure;
    for iFig = 1:nShift
        subplot(2,3,iFig);
        contourf(lowFreqRange,gammaRange,squeeze(median(squeeze(plotVar(iFig,:,:,:,:)),[1,4],'omitnan')),'lines','none');colormap jet; clim([0 1e-4])
        title([num2str(shiftLen(iFig)/1e3) ' s']); colorbar;
    end
    sgtitle(figTitle);
    % figure;
    % for iCh = 1: size(plotVar,3)
    %     subplot(size(plotVar,3),1,iCh)
    %     contourf(lowFreqRange,gammaRange,plotVar(:,:,iCh),'lines','none'); colormap jet; clim([0 1e-4])
    %     if iCh~=size(plotVar,3)
    %         xticklabels({}); yticklabels({});
    %     end
    % end
    % sgtitle(figTitle);
end

% Plot for the control/shuffled distribution
% avgShuffledModA2B = squeeze(mean(modIdxAllA2BShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
% avgShuffledModB2A = squeeze(mean(modIdxAllB2AShuffle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
% figure; contourf(lowFreqRange,gammaRange,median(avgShuffledModA2B,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;
% figure; contourf(lowFreqRange,gammaRange,median(avgShuffledModB2A,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;


% Plot for within probe
% avgModA2A = squeeze(mean(modIdxAllA2ACircle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan'));
% avgModB2B = squeeze(mean(modIdxAllB2BCircle{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
% figure; contourf(lowFreqRange,gammaRange,median(avgModA2A,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;
% figure; contourf(lowFreqRange,gammaRange,median(avgModB2B,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;

%% Organizing the data to perform statistical testing
matSize = size(modIdxAllA2B);
modIdxA2BT = reshape(modIdxAllA2B,[matSize(1)*matSize(2) 1]); modIdxA2BT(cellfun(@isempty,modIdxA2BT))=[];
modIdxB2AT = reshape(modIdxAllB2A,[matSize(1)*matSize(2) 1]); modIdxB2AT(cellfun(@isempty,modIdxB2AT))=[];
modIdxA2AT = reshape(modIdxAllA2A,[matSize(1)*matSize(2) 1]); modIdxA2AT(cellfun(@isempty,modIdxA2AT))=[];
modIdxB2BT = reshape(modIdxAllB2B,[matSize(1)*matSize(2) 1]); modIdxB2BT(cellfun(@isempty,modIdxB2BT))=[];

modIdxA2BCtrlT = reshape(modIdxAllA2BCircle,[matSize(1)*matSize(2) 1]); modIdxA2BCtrlT(cellfun(@isempty,modIdxA2BCtrlT))=[];
modIdxB2ACtrlT = reshape(modIdxAllB2ACircle,[matSize(1)*matSize(2) 1]); modIdxB2ACtrlT(cellfun(@isempty,modIdxB2ACtrlT))=[];
modIdxA2ACtrlT = reshape(modIdxAllA2ACircle,[matSize(1)*matSize(2) 1]); modIdxA2ACtrlT(cellfun(@isempty,modIdxA2ACtrlT))=[];
modIdxB2BCtrlT = reshape(modIdxAllB2BCircle,[matSize(1)*matSize(2) 1]); modIdxB2BCtrlT(cellfun(@isempty,modIdxB2BCtrlT))=[];

% Average across windows, and shuffles
nHigh    = size(gammaRange,2);
nLow     = size(lowFreqRange,2);
nSample  = size(modIdxA2BCtrlT,1);

%%
avgA2B = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxA2BT,'UniformOutput',0); avgA2B = reshape(cat(3,avgA2B{:}),[nLow*nHigh nSample]);
avgB2A = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxB2AT,'UniformOutput',0); avgB2A = reshape(cat(3,avgB2A{:}),[nLow*nHigh nSample]);
avgA2A = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxA2AT,'UniformOutput',0); avgA2A = reshape(cat(3,avgA2A{:}),[nLow*nHigh nSample]);
avgB2B = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxB2BT,'UniformOutput',0); avgB2B = reshape(cat(3,avgB2B{:}),[nLow*nHigh nSample]);

avgA2BCtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxA2BCtrlT,'UniformOutput',0); avgA2BCtrl = reshape(cat(3,avgA2BCtrl{:}),[nLow*nHigh nSample]);
avgB2ACtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxB2ACtrlT,'UniformOutput',0); avgB2ACtrl = reshape(cat(3,avgB2ACtrl{:}),[nLow*nHigh nSample]);
avgA2ACtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxA2ACtrlT,'UniformOutput',0); avgA2ACtrl = reshape(cat(3,avgA2ACtrl{:}),[nLow*nHigh nSample]);
avgB2BCtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),modIdxB2BCtrlT,'UniformOutput',0); avgB2BCtrl = reshape(cat(3,avgB2BCtrl{:}),[nLow*nHigh nSample]);

% avgA2BCtrl = cellfun(@(x)squeeze(median(x,[1 2 5],'omitnan')),modIdxA2BCtrlT,'UniformOutput',0); avgA2BCtrl = reshape(cat(3,avgA2BCtrl{:}),[nLow*nHigh nSample]);
% avgB2ACtrl = cellfun(@(x)squeeze(median(x,[1 2 5],'omitnan')),modIdxB2ACtrlT,'UniformOutput',0); avgB2ACtrl = reshape(cat(3,avgB2ACtrl{:}),[nLow*nHigh nSample]);
% avgA2ACtrl = cellfun(@(x)squeeze(median(x,[1 2 5],'omitnan')),modIdxA2ACtrlT,'UniformOutput',0); avgA2ACtrl = reshape(cat(3,avgA2ACtrl{:}),[nLow*nHigh nSample]);
% avgB2BCtrl = cellfun(@(x)squeeze(median(x,[1 2 5],'omitnan')),modIdxB2BCtrlT,'UniformOutput',0); avgB2BCtrl = reshape(cat(3,avgB2BCtrl{:}),[nLow*nHigh nSample]);


% Get the frequency pairs that are significantly different from the control
% distribution
[~,a2bPixels] = ttest(avgA2B',avgA2BCtrl'); a2bPixels = a2bPixels<0.001;
[~,b2aPixels] = ttest(avgB2A',avgB2ACtrl'); b2aPixels = b2aPixels<0.001;
[~,a2aPixels] = ttest(avgA2A',avgA2ACtrl'); a2aPixels = a2aPixels<0.001;
[~,b2bPixels] = ttest(avgB2B',avgB2BCtrl'); b2bPixels = b2bPixels<0.001;

figure;
subplot(221); imagesc(lowFreqRange,gammaRange,reshape(a2bPixels,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar; title('A--B'); axis square;
subplot(222); imagesc(lowFreqRange,gammaRange,reshape(b2aPixels,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar; title('B--A'); axis square;
subplot(223); imagesc(lowFreqRange,gammaRange,reshape(a2aPixels,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar; title('A--A'); axis square;
subplot(224); imagesc(lowFreqRange,gammaRange,reshape(b2bPixels,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar; title('B--B'); axis square;

% Combining all results
[~,allVals] = ttest([avgA2B avgB2A avgA2A avgB2B]',[avgA2BCtrl avgB2ACtrl avgA2ACtrl avgB2BCtrl]'); allVals = allVals<0.001;
[~,likeVals] = ttest([avgA2A avgB2B]',[avgA2ACtrl avgB2BCtrl]'); likeVals = likeVals<0.0003;
[~,unlikeVals] = ttest([avgA2B avgB2A]',[avgA2BCtrl avgB2ACtrl]'); unlikeVals = unlikeVals<0.0003;

figure; imagesc(lowFreqRange,gammaRange,reshape(allVals,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar;
figure; imagesc(lowFreqRange,gammaRange,reshape(likeVals,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar;
figure; imagesc(lowFreqRange,gammaRange,reshape(unlikeVals,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar;

%% Get the average PAC based on the significant pixels
avgA2B(~a2bPixels,:) = 0;
avgB2A(~b2aPixels,:) = 0;
avgA2A(~a2aPixels,:) = 0;
avgB2B(~b2bPixels,:) = 0;

figure;
subplot(221);contourf(lowFreqRange,gammaRange,reshape(avgA2B(:,5),[nHigh nLow]),'lines','none');clim([0 1e-4])  ;colorbar; colormap jet;title('A--B'); axis square;
subplot(222); contourf(lowFreqRange,gammaRange,reshape(avgB2A(:,5),[nHigh nLow]),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet; title('B--A'); axis square;
subplot(223); contourf(lowFreqRange,gammaRange,reshape(avgA2A(:,5),[nHigh nLow]),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet; title('A--A'); axis square;
subplot(224); contourf(lowFreqRange,gammaRange,reshape(avgB2B(:,5),[nHigh nLow]),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet; title('B--B'); axis square;

%% Get the maximum PAC
emptyIdx = cellfun(@isempty,(reshape(modIdxAllA2B,[matSize(1)*matSize(2) 1])));

connValsT = connValsAll; connValsT(:,1) = NaN; connValsT(:,4:end) = []; connValsT(end,:) = [];
connValsT = reshape(connValsT,[matSize(1)*matSize(2) 1]);connValsT(emptyIdx) = [];

distValT = distSitesAll; distValT(:,1) = NaN; distValT(:,4:end) = []; distValT(end,:) = [];
distValT = reshape(distValT,[matSize(1)*matSize(2) 1]); distValT(emptyIdx) = [];

avgA2B(~a2bPixels,:) = NaN;
avgB2A(~b2aPixels,:) = NaN;
avgA2A(~a2aPixels,:) = NaN;
avgB2B(~b2bPixels,:) = NaN;

%%
figure;
subplot(221); showLinearFit(connValsT,max(avgA2B,[],1,'omitnan')); box off; %xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
subplot(222); showLinearFit(connValsT,max(avgB2A,[],1,'omitnan')); box off; %xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
subplot(223); showLinearFit(connValsT,max(avgA2A,[],1,'omitnan')); box off; %xlim([-0.6 0.8]); ylim([0 5e-4]);title('A--A'); axis square;
subplot(224); showLinearFit(connValsT,max(avgB2B,[],1,'omitnan')); box off; %xlim([-0.6 0.8]); ylim([0 5e-4]);title('B--B'); axis square;

%% Divide into frequency pairs
avgA2BT = reshape(avgA2B,[nHigh nLow nSample]);
avgB2AT = reshape(avgB2A,[nHigh nLow nSample]);
avgA2AT = reshape(avgA2A,[nHigh nLow nSample]);
avgB2BT = reshape(avgB2B,[nHigh nLow nSample]);

for iPlot = 1:4
    switch iPlot
        case 1
            a2b = median(avgA2BT(1:6,1:2,:),[1 2],'omitnan');
            b2a = median(avgB2AT(1:6,1:2,:),[1 2],'omitnan');
            a2a = median(avgA2AT(1:6,1:2,:),[1 2],'omitnan');
            b2b = median(avgB2BT(1:6,1:2,:),[1 2],'omitnan');
            pltTitle = 'Theta to Low Gamma';
        case 2
            a2b = median(avgA2BT(7:11,1:2,:),[1 2],'omitnan');
            b2a = median(avgB2AT(7:11,1:2,:),[1 2],'omitnan');
            a2a = median(avgA2AT(7:11,1:2,:),[1 2],'omitnan');
            b2b = median(avgB2BT(7:11,1:2,:),[1 2],'omitnan');
            pltTitle = 'Theta to High Gamma';
        case 3
            a2b = median(avgA2BT(1:6,3:4,:),[1 2],'omitnan');
            b2a = median(avgB2AT(1:6,3:4,:),[1 2],'omitnan');
            a2a = median(avgA2AT(1:6,3:4,:),[1 2],'omitnan');
            b2b = median(avgB2BT(1:6,3:4,:),[1 2],'omitnan');
            pltTitle = 'Alpha to Low Gamma';
        case 4
            a2b = mean(avgA2BT(7:11,3:4,:),[1 2],'omitnan');
            b2a = mean(avgB2AT(7:11,3:4,:),[1 2],'omitnan');
            a2a = mean(avgA2AT(7:11,3:4,:),[1 2],'omitnan');
            b2b = mean(avgB2BT(7:11,3:4,:),[1 2],'omitnan');
            pltTitle = 'Alpha to High Gamma';
    end
    % figure;
    % subplot(221); showLinearFit(connValsT,squeeze(a2b),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
    % subplot(222); showLinearFit(connValsT,squeeze(b2a),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
    % subplot(223); showLinearFit(connValsT,squeeze(a2a),0.3,1.85e-4,1.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('A--A'); axis square;
    % subplot(224); showLinearFit(connValsT,squeeze(a2b),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('B--B'); axis square;
    % xlabel('Functional connectivity'); ylabel('Modulation index')
    % sgtitle(pltTitle);

    figure;
    boxplot([squeeze(a2b) squeeze(b2a) squeeze(a2a) squeeze(b2b)],{'A-B','B-A','A-A','B-B'}); ylim([0 2e-4]);
    title(pltTitle); box off;
end
%%  Laminar analysis
% Dividing and averaging within compartments
for iVar = 1:4
    switch iVar
        case 1
            var = modIdxA2BT;
        case 2
            var = modIdxB2AT;
        case 3
            var = modIdxA2AT;
        case 4
            var = modIdxB2BT;
    end
    clear super mid deep
    super = cellfun(@(x) reshape(squeeze(median(x(:,:,:,1:6),[1 4])),[nLow*nHigh 1]),var,'UniformOutput',0);    super = cat(2,super{:});
    mid   = cellfun(@(x) reshape(squeeze(median(x(:,:,:,7:12),[1 4])),[nLow*nHigh 1]),var,'UniformOutput',0);   mid = cat(2,mid{:});
    deep  = cellfun(@(x) reshape(squeeze(median(x(:,:,:,13:end),[1 4])),[nLow*nHigh 1]),var,'UniformOutput',0); deep = cat(2,deep{:});

    switch iVar
        case 1
            a2bLaminar = cat(3,super,mid,deep); a2bLaminar = permute(a2bLaminar,[1 3 2]);
            a2bLaminar(~a2bPixels,:,:) = NaN; a2bLaminar = reshape(a2bLaminar,[nHigh nLow 3 nSample]);
        case 2
            b2aLaminar = cat(3,super,mid,deep); b2aLaminar = permute(b2aLaminar,[1 3 2]);
            b2aLaminar(~a2bPixels,:,:) = NaN; b2aLaminar = reshape(b2aLaminar,[nHigh nLow 3 nSample]);
        case 3
            a2aLaminar = cat(3,super,mid,deep); a2aLaminar = permute(a2aLaminar,[1 3 2]);
            a2aLaminar(~a2bPixels,:,:) = NaN; a2aLaminar = reshape(a2aLaminar,[nHigh nLow 3 nSample]);
        case 4
            b2bLaminar = cat(3,super,mid,deep); b2bLaminar = permute(b2bLaminar,[1 3 2]);
            b2bLaminar(~a2bPixels,:,:) = NaN; b2bLaminar = reshape(b2bLaminar,[nHigh nLow 3 nSample]);
    end
end

%% Plot the compartments
for iPlot = 1%:4
    switch iPlot
        case 1
            a2b = squeeze(median(a2bLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            b2a = squeeze(median(b2aLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            a2a = squeeze(median(a2aLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            b2b = squeeze(median(b2bLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            pltTitle = 'Theta to Low Gamma';
        case 2
            a2b = squeeze(median(a2bLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            b2a = squeeze(median(b2aLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            a2a = squeeze(median(a2aLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            b2b = squeeze(median(b2bLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            pltTitle = 'Theta to High Gamma';
        case 3
            a2b = squeeze(median(a2bLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            b2a = squeeze(median(b2aLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            a2a = squeeze(median(a2aLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            b2b = squeeze(median(b2bLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            pltTitle = 'Alpha to Low Gamma';
        case 4
            a2b = squeeze(median(a2bLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            b2a = squeeze(median(b2aLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            a2a = squeeze(median(a2aLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            b2b = squeeze(median(b2bLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            pltTitle = 'Alpha to High Gamma';
    end

    for iLayer = 1:3
        switch iLayer
            case 1
                layerTitle = 'Superficial';
            case 2
                layerTitle = 'Middle';
            case 3
                layerTitle = 'Deep';
        end

        % figure;
        % subplot(221); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
        % subplot(222); showLinearFit(connValsT,squeeze(b2a(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
        % subplot(223); showLinearFit(connValsT,squeeze(a2a(iLayer,:)),0.3,1.85e-4,1.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('A--A'); axis square;
        % subplot(224); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('B--B'); axis square;
        % xlabel('Functional connectivity'); ylabel('Modulation index')
        % sgtitle([pltTitle '- ' layerTitle]);

        figure;
        boxplot([squeeze(a2b(iLayer,:)) ;squeeze(b2a(iLayer,:)); squeeze(a2a(iLayer,:)); squeeze(b2b(iLayer,:))]',{'A-B','B-A','A-A','B-B'});
        sgtitle([pltTitle '- ' layerTitle]); box off; ylim([0 2e-4]);
    end
end

%% Combining A-B with B-A and A-A with B-B
clear modIdxUnlikePairs modIdxLikePairs modIdxUnlikePairsCtrl modIdxLikePairsCtrl
% Combining the like (within probe) and unlike (between probe) comparisons
modIdxUnlikePairs = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')), cellfun(@(x,y)(x+y)./2, modIdxA2BT, modIdxB2AT,'un',0),'un',0);
modIdxUnlikePairs = reshape(cat(3,modIdxUnlikePairs{:}),[nLow*nHigh nSample]);

modIdxLikePairs = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),[modIdxA2AT ;modIdxB2BT],'un',0);
modIdxLikePairs = reshape(cat(3,modIdxLikePairs{:}),[nLow*nHigh nSample*2]);

modIdxUnlikePairsCtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),cellfun(@(x,y)(x+y)./2, modIdxA2BCtrlT, modIdxB2ACtrlT,'un',0),'un',0);
modIdxUnlikePairsCtrl = reshape(cat(3,modIdxUnlikePairsCtrl{:}),[nLow*nHigh nSample]);

modIdxLikePairsCtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),[modIdxA2ACtrlT ;modIdxB2BCtrlT],'un',0);
modIdxLikePairsCtrl = reshape(cat(3,modIdxLikePairsCtrl{:}),[nLow*nHigh nSample*2]);

[~,unlikePixels] = ttest(modIdxUnlikePairs',modIdxUnlikePairsCtrl'); unlikePixels = unlikePixels<0.0003;
[~,likePixels]   = ttest(modIdxLikePairs',modIdxLikePairsCtrl'); likePixels = likePixels<0.0003;
figure;
subplot(121);imagesc(lowFreqRange,gammaRange,reshape(unlikePixels,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar; axis square;
subplot(122); imagesc(lowFreqRange,gammaRange,reshape(likePixels,[nHigh nLow]));set(gca,'YDir','normal');colormap gray; colorbar; axis square;

% Get the average MI based on the significant pixels
modIdxLikePairs(~likePixels,:)     = 0;
modIdxUnlikePairs(~unlikePixels,:) = 0;

figure;
subplot(121); contourf(lowFreqRange,gammaRange,reshape(modIdxUnlikePairs(:,4),[nHigh nLow]),'lines','none');clim([0 1e-4])  ;colorbar; colormap jet;title('Between probes'); axis square;
subplot(122); contourf(lowFreqRange,gammaRange,reshape(modIdxLikePairs(:,4),[nHigh nLow]),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet; title('Within probe'); axis square;

%%
modIdxLikePairs = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),[modIdxA2AT ;modIdxB2BT],'un',0);
modIdxLikePairs = reshape(cat(3,modIdxLikePairs{:}),[nLow*nHigh nSample*2]);

modIdxUnlikePairsCtrl = cellfun(@(x)squeeze(median(x,[1 4],'omitnan')),cellfun(@(x,y)(x+y)./2, modIdxA2BCtrlT, modIdxB2ACtrlT,'un',0),'un',0);
modIdxUnlikePairsCtrl = reshape(cat(3,modIdxUnlikePairsCtrl{:}),[nLow*nHigh nSample]);

modIdxLikePairs(~likePixels,:)         = NaN;
modIdxUnlikePairs(~unlikePixels,:)     = NaN;
modIdxLikePairsCtrl(~likePixels,:)     = NaN;
modIdxUnlikePairsCtrl(~unlikePixels,:) = NaN;

modIdxLikePairsT     = reshape(modIdxLikePairs,[nHigh nLow nSample*2]);
modIdxLikePairsCtrlT = reshape(modIdxLikePairsCtrl,[nHigh nLow nSample*2]);

modIdxUnlikePairsT     = reshape(modIdxUnlikePairs,[nHigh nLow nSample]);
modIdxUnlikePairsCtrlT = reshape(modIdxUnlikePairsCtrl,[nHigh nLow nSample]);

%%
[yIndex,edgeVals] = discretize(connValsT,-0.4:0.2:0.8);
loc1 = yIndex==1 | yIndex==2 | yIndex==3; loc2 = yIndex==4 | yIndex==5 | yIndex==6;

for iPlot = 1:4
    switch iPlot
        case 1
            unlikePairs     = squeeze(median(modIdxUnlikePairsT(1:6,1:2,:),[1 2],'omitnan'));
            unlikePairsCtrl = squeeze(median(modIdxUnlikePairsCtrlT(1:6,1:2,:),[1 2],'omitnan'));

            likePairs     = squeeze(median(modIdxLikePairsT(1:6,1:2,:),[1 2],'omitnan'));
            likePairsCtrl = squeeze(median(modIdxLikePairsCtrlT(1:6,1:2,:),[1 2],'omitnan'));

            pltTitle = 'Theta to Low Gamma';
        case 2
            unlikePairs     = squeeze(median(modIdxUnlikePairsT(7:11,1:2,:),[1 2],'omitnan'));
            unlikePairsCtrl = squeeze(median(modIdxUnlikePairsCtrlT(7:11,1:2,:),[1 2],'omitnan'));

            likePairs     = squeeze(median(modIdxLikePairsT(7:11,1:2,:),[1 2],'omitnan'));
            likePairsCtrl = squeeze(median(modIdxLikePairsCtrlT(7:11,1:2,:),[1 2],'omitnan'));

            pltTitle = 'Theta to High Gamma';
        case 3
            unlikePairs     = squeeze(median(modIdxUnlikePairsT(1:6,3:4,:),[1 2],'omitnan'));
            unlikePairsCtrl = squeeze(median(modIdxUnlikePairsCtrlT(1:6,3:4,:),[1 2],'omitnan'));

            likePairs     = squeeze(median(modIdxLikePairsT(1:6,3:4,:),[1 2],'omitnan'));
            likePairsCtrl = squeeze(median(modIdxLikePairsCtrlT(1:6,3:4,:),[1 2],'omitnan'));

            pltTitle = 'Alpha to Low Gamma';

        case 4
            unlikePairs     = squeeze(median(modIdxUnlikePairsT(7:11,3:4,:),[1 2],'omitnan'));
            unlikePairsCtrl = squeeze(median(modIdxUnlikePairsCtrlT(7:11,3:4,:),[1 2],'omitnan'));

            likePairs     = squeeze(median(modIdxLikePairsT(7:11,3:4,:),[1 2],'omitnan'));
            likePairsCtrl = squeeze(median(modIdxLikePairsCtrlT(7:11,3:4,:),[1 2],'omitnan'));

            pltTitle = 'Alpha to High Gamma';
    end
    % figure;
    % subplot(221); showLinearFit(connValsT,squeeze(a2b),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
    % subplot(222); showLinearFit(connValsT,squeeze(b2a),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
    % subplot(223); showLinearFit(connValsT,squeeze(a2a),0.3,1.85e-4,1.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('A--A'); axis square;
    % subplot(224); showLinearFit(connValsT,squeeze(a2b),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('B--B'); axis square;
    % xlabel('Functional connectivity'); ylabel('Modulation index')
    % sgtitle(pltTitle);

    % figure;
    % boxplot([[unlikePairs; NaN(25,1)] likePairs [unlikePairsCtrl; NaN(25,1)] likePairsCtrl],{'Unlike pairs','Like pairs','Control-unlike','control-like'}); ylim([0 2e-4]);
    figure;
    histogram(unlikePairs(loc1),0:0.05e-4:1e-4); hold on;
    histogram(unlikePairs(loc2),0:0.05e-4:1e-4); hold on;
    xlim([0 1e-4]); ylim([0 15]);box off;
    legend('FC<0.2','FC>0.2','Location','northeast');

    % figure;
    % histogram(likePairs([loc1;loc1]),0:0.05e-4:1e-4); hold on;
    % histogram(likePairs([loc2; loc2]),0:0.05e-4:1e-4); hold on;
    % xlim([0 1e-4]); ylim([0 8]);box off;
    % legend('FC<0.2','FC>0.2','Location','northeast');

    title(pltTitle); box off;
end

%% Laminar
% Dividing and averaging within compartments
for iVar = 1:2
    switch iVar
        case 1
            var =  cellfun(@(x,y)(x+y)./2, modIdxA2BT, modIdxB2AT,'un',0); % Unlike pairs
        case 2
            var = [modIdxA2AT ;modIdxB2BT];
        case 3
            var = cellfun(@(x,y)(x+y)./2, modIdxA2BCtrlT, modIdxB2ACtrlT,'un',0);
        case 4
            var = [modIdxA2ACtrlT; modIdxB2VtrlT];
    end
    clear super mid deep
    super = cellfun(@(x) reshape(squeeze(median(x(:,:,:,1:6),[1 4])),[nLow*nHigh 1]),var,'UniformOutput',0);    super = cat(2,super{:});
    mid   = cellfun(@(x) reshape(squeeze(median(x(:,:,:,7:12),[1 4])),[nLow*nHigh 1]),var,'UniformOutput',0);   mid = cat(2,mid{:});
    deep  = cellfun(@(x) reshape(squeeze(median(x(:,:,:,13:end),[1 4])),[nLow*nHigh 1]),var,'UniformOutput',0); deep = cat(2,deep{:});

    switch iVar
        case 1
            unlikeLaminar = cat(3,super,mid,deep); unlikeLaminar = permute(unlikeLaminar,[1 3 2]);
            unlikeLaminar(~unlikePixels,:,:) = NaN; unlikeLaminar = reshape(unlikeLaminar,[nHigh nLow 3 nSample]);
        case 2
            likeLaminar = cat(3,super,mid,deep); likeLaminar = permute(likeLaminar,[1 3 2]);
            likeLaminar(~likePixels,:,:) = NaN; likeLaminar = reshape(likeLaminar,[nHigh nLow 3 nSample*2]);
        case 3
            unlikeCtrlLaminar = cat(3,super,mid,deep); unlikeCtrlLaminar = permute(unlikeCtrlLaminar,[1 3 2]);
            unlikeCtrlLaminar(~unlikePixels,:,:) = NaN; unlikeCtrlLaminar = reshape(unlikeCtrlLaminar,[nHigh nLow 3 nSample]);
        case 4
            likeCtrlLaminar = cat(3,super,mid,deep); likeCtrlLaminar = permute(likeCtrlLaminar,[1 3 2]);
            likeCtrlLaminar(~likePixels,:,:) = NaN; likeCtrlLaminar = reshape(likeCtrlLaminar,[nHigh nLow 3 nSample*2]);
    end
end

%% Plot the compartments
for iPlot = 1:4
    switch iPlot
        case 1
            unlikeLayers = squeeze(median(unlikeLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            likeLayers   = squeeze(median(likeLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            % a2a = squeeze(median(a2aLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            % b2b = squeeze(median(b2bLaminar(1:6,1:2,:,:),[1 2],'omitnan'));
            pltTitle = 'Theta to Low Gamma';
        case 2
            unlikeLayers = squeeze(median(unlikeLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            likeLayers   = squeeze(median(likeLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            % a2a = squeeze(median(a2aLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            % b2b = squeeze(median(b2bLaminar(7:11,1:2,:,:),[1 2],'omitnan'));
            pltTitle = 'Theta to High Gamma';
        case 3
            unlikeLayers = squeeze(median(unlikeLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            likeLayers   = squeeze(median(likeLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            % a2a = squeeze(median(a2aLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            % b2b = squeeze(median(b2bLaminar(1:6,3:4,:,:),[1 2],'omitnan'));
            pltTitle = 'Alpha to Low Gamma';
        case 4
            unlikeLayers = squeeze(median(unlikeLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            likeLayers   = squeeze(median(likeLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            % a2a = squeeze(median(a2aLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            % b2b = squeeze(median(b2bLaminar(7:11,3:4,:,:),[1 2],'omitnan'));
            pltTitle = 'Alpha to High Gamma';
    end

    figure;
    for iP = 1:3
        switch iP
            case 1
                layerTitle = 'Superficial';
            case 2
                layerTitle = 'Middle';
            case 3
                layerTitle = 'Deep';
        end
        subplot(1,3,iP);
        histogram(unlikeLayers(iP,loc1),0:0.05e-4:1e-4); hold on;
        histogram(unlikeLayers(iP,loc2),0:0.05e-4:1e-4);
        xlim([0 1e-4]); ylim([0 8]);box off;
        legend('FC<0.2','FC>0.2','Location','northeast');
        title(layerTitle);axis square;
    end
    % histogram(unlikePairs(iP,loc2),0:0.05e-4:1e-4); hold on;
    % xlim([0 1e-4]); ylim([0 8]);box off;
    % legend('FC<0.2','FC>0.2','Location','northeast');

    sgtitle(pltTitle); box off;

    % figure;
    % subplot(121); boxplot(unlikeLayers',{'Superficial','Middle','Deep'}); ylim([0 2e-4]); title('Unlike pairs'); box off; axis square;
    % subplot(122); boxplot(likeLayers',{'Superficial','Middle','Deep'});ylim([0 2e-4]);  title('Like pairs'); box off; axis square;
    % sgtitle(pltTitle);
    % figure;
    % for iLayer = 1:3
    %     switch iLayer
    %         case 1
    %             layerTitle = 'Superficial';
    %         case 2
    %             layerTitle = 'Middle';
    %         case 3
    %             layerTitle = 'Deep';
    %     end
    %
    %     % figure;
    %     % subplot(221); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('A--B'); axis square;
    %     % subplot(222); showLinearFit(connValsT,squeeze(b2a(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 1e-4]);title('B--A'); axis square;
    %     % subplot(223); showLinearFit(connValsT,squeeze(a2a(iLayer,:)),0.3,1.85e-4,1.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('A--A'); axis square;
    %     % subplot(224); showLinearFit(connValsT,squeeze(a2b(iLayer,:)),0.3,0.85e-4,0.8e-4); box off; xlim([-0.6 0.8]); ylim([0 2e-4]);title('B--B'); axis square;
    %     % xlabel('Functional connectivity'); ylabel('Modulation index')
    %     % sgtitle([pltTitle '- ' layerTitle]);
    %
    %     subplot(1,3,iLayer);
    %     boxplot([[squeeze(unlikeLayers(iLayer,:)) NaN(1,25)] ;squeeze(likeLayers(iLayer,:))]',{'Unlike pairs','Like pairs'});
    %     title(layerTitle); box off; ylim([0 2e-4]); axis square;
    % end
    % sgtitle(pltTitle);
end



%% Dividing the data in 10 s increments
dataLen = size(probeA,1);
stepLen = 30e3;
winLen  = 30e3:stepLen:dataLen;

clear modIdx
lowFreqB  = single(eegfilt(probeB(:,chB(1):chB(2))',fs,thetaBand(1),thetaBand(2)))';
lowFreqA  = single(eegfilt(probeA(:,chA(1):chA(2))',fs,thetaBand(1),thetaBand(2)))';
highFreqB = single(eegfilt(probeB(:,chB(1):chB(2))',fs,gammaBand(1),gammaBand(2)))';
highFreqA = single(eegfilt(probeA(:,chA(1):chA(2))',fs,gammaBand(1),gammaBand(2)))';

for iWin = 1:size(winLen,2)
    % Divide the data into chunks
    winSize = 1:winLen(iWin);
    iCh = 1;
    while ~(winSize(1)>dataLen)
        if winSize(end)>dataLen; winSize = winSize(1):dataLen;end
        if ~(dataLen-winSize(1)<10e3)
            [modIdxA2B(iCh,iWin,:),~] = getPhaseAmpCoupling(lowFreqB(winSize,:),highFreqA(winSize,:),0);
            [modIdxB2A(iCh,iWin,:),~] = getPhaseAmpCoupling(lowFreqA(winSize,:),highFreqB(winSize,:),0);
            [modIdxAA(iCh,iWin,:),~]  = getPhaseAmpCoupling(lowFreqA(winSize,:),highFreqA(winSize,:),0);
            [modIdxBB(iCh,iWin,:),~]  = getPhaseAmpCoupling(lowFreqB(winSize,:),highFreqB(winSize,:),0);
        end
        iCh = iCh+1;
        winSize = winSize+(winLen(iWin)/2);
    end
end

for iMod = 1:4
    clear modIdx meanModIdx
    switch iMod
        case 1
            modIdx = modIdxA2B;
            figTitle = 'A-B';
        case 2
            modIdx = modIdxB2A;
            figTitle = 'B-A';
        case 3
            modIdx = modIdxAA;
            figTitle = 'A-A';
        case 4
            modIdx = modIdxBB;
            figTitle = 'B-B';
    end

    modIdx(modIdx==0) = NaN;
    meanModIdx = squeeze(median(modIdx,1,'omitnan'));
    stdModIdx  = squeeze(mad(modIdx,1,1));
    figure;
    for iP = 1:size(meanModIdx,2)
        subplot(5,5,iP);
        if iP == 1; ylabel('Modulation index'); xlabel('Data length (s)'); end
        plot(winLen,meanModIdx(:,iP),'k','LineWidth',1.5); hold on;
        plot(winLen,meanModIdx(:,iP)-2.*stdModIdx(:,iP),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
        plot(winLen,meanModIdx(:,iP)+2.*stdModIdx(:,iP),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
        xticks(0:100e3:dataLen); xticklabels(0:100:dataLen/1e3); box off; ylim([0 2e-3]);
        yline(0)
    end
    sgtitle(figTitle);
end

%% Functions
function [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
    probeA,probeB,badElecA,badElecB,badTimes,chA,chB,gammaRange,lowFreqRange)

fs = 1e3;
if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
    varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
    varIdx  = ismember('amplitudeA',varInfo); % Check if the variable exists
else
    % Check if the mat file exists already
    if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
        varIdx = 0;
    else
        varIdx = 1;
    end
end


if ~varIdx
    % Remove bad channelscommonDir
    probeA(:,badElecA) = [];
    probeB(:,badElecB) = [];

    % Remove bad times
    probeA(badTimes,:) = [];
    probeB(badTimes,:) = [];

    dataLen = fix(size(probeA,1)/1000)*1000;

    % if chA(1)== 0 || chB(1)==0; break; end

    % Get the bandlimited data
    clear modA2B modB2A

    % High frequencies
    disp('Getting high frequency amplitudes ....');

    for iHigh = 1:size(gammaRange,2)
        fLow = gammaRange(iHigh);
        fHigh = gammaRange(iHigh)+5; % 5 Hz bandwidth

        highFreqA = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh))';
        highFreqB = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh))';

        amplitudeA(iHigh,:,:) = abs(hilbert(highFreqA));
        amplitudeB(iHigh,:,:) = abs(hilbert(highFreqB));

    end


    % Low frequencies
    disp('Filtering low frequency signals....');
    for iLow = 1:size(lowFreqRange,2)
        fLow = lowFreqRange(iLow);
        fHigh = lowFreqRange(iLow)+2; % 2 Hz bandwidth

        lowFreqA  = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh))';
        lowFreqB = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh))';

        phaseA(iLow,:,:)   = angle(hilbert(lowFreqA)); % Get phase of low frequency
        phaseB(iLow,:,:)   = angle(hilbert(lowFreqB)); % Get phase of low frequency

    end

    if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
        save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
            'amplitudeA','amplitudeB','phaseA','phaseB');

    else
        save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
            'amplitudeA','amplitudeB','phaseA','phaseB','-append');

    end

else
    vars = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
    amplitudeA = vars.amplitudeA;
    amplitudeB = vars.amplitudeB;
    phaseA     = vars.phaseA;
    phaseB     = vars.phaseB;
end
end
%%
function [phaseAShuff, phaseBShuff] = calculateRandPhaseSignal(monkeyName,expDate,hemisphere,fileNum,...
    probeA,probeB,badElecA,badElecB,badTimes,chA,chB,lowFreqRange)

fs = 1e3;
% if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
%     varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
%     varIdx  = ismember('amplitudeA',varInfo); % Check if the variable exists
% else
%     % Check if the mat file exists already
%     if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
%         varIdx = 0;
%     else
%         varIdx = 1;
%     end
% end

% if ~varIdx
    % Remove bad channels
    probeA(:,badElecA) = [];
    probeB(:,badElecB) = [];

    % Remove bad times
    probeA(badTimes,:) = [];
    probeB(badTimes,:) = [];

    dataLen = fix(size(probeA,1)/1000)*1000;

    % if chA(1)== 0 || chB(1)==0; break; end

    % Get the bandlimited data
    clear modA2B modB2A

    % % High frequencies
    % disp('Getting high frequency amplitudes ....');
    % 
    % for iHigh = 1:size(gammaRange,2)
    %     fLow = gammaRange(iHigh);
    %     fHigh = gammaRange(iHigh)+5; % 5 Hz bandwidth
    % 
    %     highFreqA = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh))';
    %     highFreqB = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh))';
    % 
    %     amplitudeA(iHigh,:,:) = abs(hilbert(highFreqA));
    %     amplitudeB(iHigh,:,:) = abs(hilbert(highFreqB));
    % 
    % end


    % Low frequencies
    disp('Filtering low frequency signals....');
    for iLow = 1:size(lowFreqRange,2)
        fLow = lowFreqRange(iLow);
        fHigh = lowFreqRange(iLow)+2; % 2 Hz bandwidth

        lowFreqA  = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh))';
        lowFreqB = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh))';

        % for iRep = 1:5
            rng('shuffle');
            comb1 = randperm(dataLen);

            fftA   = fft(lowFreqA); % Fourier transform
            magA   = abs(fftA); % Magnitude
            phaseA = angle(fftA); % Phase
            newPhaseA = phaseA(comb1,:);
            newLowA    = real(ifft(magA(1:size(newPhaseA,1),:).*newPhaseA));

            fftB   = fft(lowFreqB);
            magB   = abs(fftB);
            phaseB = angle(fftB);
            newPhaseB = phaseB(comb1,:);
            newLowB    = real(ifft(magB(1:size(newPhaseB,1),:).*newPhaseB));

            phaseAShuff(iLow,:,:)   = angle(hilbert(newLowA)); % Get phase of low frequency
            phaseBShuff(iLow,:,:)   = angle(hilbert(newLowB)); % Get phase of low frequency
        % end

    end

%     if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
%         save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
%             'amplitudeA','amplitudeB','phaseA','phaseB');
% 
%     else
%         save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
%             'amplitudeA','amplitudeB','phaseA','phaseB','-append');
% 
%     end
% 
% else
%     vars = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
%     amplitudeA = vars.amplitudeA;
%     amplitudeB = vars.amplitudeB;
%     phaseA     = vars.phaseA;
%     phaseB     = vars.phaseB;
% end
end

%% Function to fit a line
function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
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
text(textLocX, textLocY2,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]);
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

if nargin<3
    textLocX  = max(xVal)-0.2*max(xVal);
    textLocY1 = max(yVal)-0.2*max(yVal);
    textLocY2 = max(yVal)-0.3*max(yVal);
end

text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end