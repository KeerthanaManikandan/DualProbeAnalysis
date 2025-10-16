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

iM = 2; % 1 - Charlie Sheen, 2 - Whiskey
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
[allVars] = getDualProbeCorrelations(monkeyName, hemisphere, allDates, datFileNumAll,allProbeData,...
    badElecB,estChInCortexA,estChInCortexB,connValsAll,distSitesAll,goodRuns);

% Check plotting...


%% Perform Phase amplitude coupling 
% Get the PAC comodulogram
clear modIdxAllA2B

thetaVals  = thetaBand;

alphaRange = alphaBand(1):2:alphaBand(2);
% alphaRange = alphaRange(2:end);
% alphaVals  = [alphaRange-1; alphaRange+1]';

betaRange = betaBand(1):2:betaBand(2);
% betaRange = betaRange(2:end);
% betaVals  = [betaRange-1; betaRange+1]';

gammaRange = gammaBand(1):5:gammaBand(2);
gammaRange(gammaRange>=60 & gammaRange<=65)=[];
% gammaRange = gammaRange(2:end);
% gammaVals = [gammaRange-2; gammaRange+2]';

lowFreqRange = 6:2:30; 

for iDate = 2:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)        
        fileNum = datFileNum(iRun);
        clear probeA probeB chA chB

        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        % Get the ephys data
        clear probeA probeB rawA rawB
        probeA = allProbeData{fileNum,iDate}.probe1Ch;
        probeB = allProbeData{fileNum,iDate}.probe2Ch;

        % Remove bad channels
        probeA(:,badElecA{fileNum,iDate}) = [];
        probeB(:,badElecB{fileNum,iDate}) = [];

        % Remove bad times
        probeA(allBadTimes{fileNum,iDate},:) = [];
        probeB(allBadTimes{fileNum,iDate},:) = [];

        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);
        dataLen = fix(size(probeA,1)/1000)*1000;

        if chA(1)== 0 || chB(1)==0; continue; end

        for iHigh = 1:size(gammaRange,2)
            fLow = gammaRange(iHigh);
            fHigh = gammaRange(iHigh)+10;
            highFreqA(iHigh,:,:) = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,1e3))';
            highFreqB(iHigh,:,:) = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,1e3))';
        end

        for iLow = 1:size(lowFreqRange,2)
            fLow = lowFreqRange(iLow);
            fHigh = lowFreqRange(iLow)+10;
            filtOrder = 3*fs/lowFreqRange(iLow);
            epochLen = round(3*filtOrder);
            
            if ~mod(dataLen,epochLen)==0
                if ~isprime(dataLen/1e3)
                    epochLen = ceil(epochLen/1000)*1e3;
                else
                    epochLen=0;
                end
            end
            lowFreqA(iLow,:,:) = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,epochLen))';
            lowFreqB(iLow,:,:) = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,epochLen))';
        end

        % Get comodulogram
        for iHigh = 1:size(gammaRange,2)
            iWin = 1;
            winLen  = 1:200e3; % Window size to calculate modulation index
            stepLen = 100e3; % Step size

            while ~(winLen(end)>dataLen)
                for iLow = 1:size(lowFreqRange,2)
                    [modIdxAllA2B{iRun,iDate}(iWin,iHigh,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowFreqB(iLow,winLen,:)),squeeze(highFreqA(iHigh,winLen,:)),0);
                    [modIdxAllB2A{iRun,iDate}(iWin,iHigh,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowFreqA(iLow,winLen,:)),squeeze(highFreqB(iHigh,winLen,:)),0);
                end
                winLen = winLen+stepLen;
                iWin = iWin+1;
            end
        end



        % for iHigh = 1:size(gammaVals,1)
        %     % Filter High frequencies
        %     clear highFreqA highFreqB
        %     highFreqA = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,gammaVals(iHigh,1),gammaVals(iHigh,2),1e3))';
        %     highFreqB = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,gammaVals(iHigh,1),gammaVals(iHigh,2),1e3))';
        % 
        %     for iFreq = 1:3
        %         clear freqVal
        %         switch iFreq
        %             case 1
        %                 freqVal = thetaVals;
        %             case 2
        %                 freqVal = alphaVals;
        %             case 3
        %                 freqVal = betaVals;
        %         end
        % 
        %         disp(['High: ' num2str(iHigh) ' Freq: ' num2str(iFreq)]);
        % 
        %         for iLow = 1:size(freqVal,1)
        %             clear lowFreqB lowFreqA
        %             filtOrder = 3*fs/freqVal(iLow,1);
        %             epochLen = round(3*filtOrder);
        % 
        %             if ~mod(dataLen,epochLen)==0
        %                 if ~isprime(dataLen/1e3)
        %                     epochLen = ceil(epochLen/1000)*1e3;
        %                 else
        %                     epochLen=0;
        %                 end
        %             end
        % 
        %             lowFreqA = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,freqVal(iLow,1),freqVal(iLow,2),epochLen))';
        %             lowFreqB = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,freqVal(iLow,1),freqVal(iLow,2),epochLen))';
        % 
        %             % Split the data into 200 s intervals, and then
        %             % calculate Modulation index
        %             iWin = 1;
        %             winLen  = 1:200e3; % Window size to calculate modulation index
        %             stepLen = 100e3; % Step size
        % 
        %             while ~(winLen(end)>dataLen)
        %                 [modIdxAllA2B{iRun,iDate}(iFreq,iWin,iHigh,iLow,:),~] = getPhaseAmpCoupling(lowFreqB(winLen,:),highFreqA(winLen,:),0);
        %                 [modIdxAllB2A{iRun,iDate}(iFreq,iWin,iHigh,iLow,:),~] = getPhaseAmpCoupling(lowFreqA(winLen,:),highFreqB(winLen,:),0);
        % 
        %                 % Shuffle the amplitude and calculate the
        %                 % modulation index -- Null distribution
        %                 % for iShuffle = 1:10
        %                 %     rng('shuffle');                           
        %                 %     highFreqShuffle =  highFreqA(randperm(size(highFreqA,1)),:);
        %                 %     [modIdxShuffle(iFreq,iWin,iShuffle,iHigh,iLow,:),~] = getPhaseAmpCoupling(lowFreqB(winLen,:),highFreqShuffle(winLen,:),0);
        %                 % end
        % 
        %                 winLen = winLen+stepLen;
        %                 iWin = iWin+1;
        %             end   
        %         end
        %     end
        % end

        % % Rearranging MI
        % thetaCols = squeeze(modIdxAllA2B(1,:,:,:));
        % modulogramVals = cat(2,thetaCols(:,:,1:size(thetaVals,1),:),squeeze(modIdxAllA2B(2,:,:,1:size(alphaVals,1),:)), ...
        %     squeeze(modIdxAllA2B(3,:,:,1:size(betaVals,1),:)));

    end
end

%% Dividing the data in 10 s increments
dataLen = size(probeA,1);
stepLen = 10e3; 
winLen  = 10e3:stepLen:dataLen;
clear modIdx
for iWin = 1:size(winLen,2)
    % Divide the data into chunks
    winSize = 1: winLen(iWin);
    iCh = 1;
    while ~(winSize(1)>dataLen)
        if winSize(end)>dataLen; winSize = winSize(1):dataLen;end
        if ~(dataLen-winSize(1)<10e3)
            lowFreqB = pATheta(winSize,:);
            highFreqA = pAGamma(winSize,:);
            [modIdx(iCh,iWin,:),~] = getPhaseAmpCoupling(lowFreqB,highFreqA,0);
        end
        iCh = iCh+1;
        winSize = winSize+(winLen(iWin)/2);
    end
end

modIdx(modIdx==0) = NaN;
meanModIdx = squeeze(median(modIdx,1,'omitnan'));
stdModIdx  = squeeze(mad(modIdx,1,1));
figure;
for iP = 1:size(meanModIdx,2)
    subplot(5,5,iP); 
    plot(winLen,meanModIdx(:,iP),'k'); hold on;
    plot(winLen,meanModIdx(:,iP)-2.*stdModIdx(:,iP),'--','Color',[0.5 0.5 0.5]);
    plot(winLen,meanModIdx(:,iP)+2.*stdModIdx(:,iP),'--','Color',[0.5 0.5 0.5]);
    xticks(0:100e3:dataLen); xticklabels(0:100:dataLen/1e3); box off; ylim([0 2e-3]);
    yline(0);
end

%% Comodulogram
