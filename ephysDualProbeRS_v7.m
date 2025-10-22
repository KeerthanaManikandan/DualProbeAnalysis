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
gammaRange = gammaBand(1):5:gammaBand(2);
gammaRange(gammaRange>=60 & gammaRange<=65)=[];

lowFreqRange = 6:2:30; 

for iDate = 2:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);

        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
            clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

            chA = estChInCortexA{iDate}(iRun,:);
            chB = estChInCortexB{iDate}(iRun,:);

            if chA(1)== 0 || chB(1)==0; continue; end

            [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
                allProbeData{fileNum,iDate}.probe1Ch,allProbeData{fileNum,iDate}.probe2Ch,badElecA{fileNum,iDate},...
                badElecB{fileNum,iDate},allBadTimes{fileNum,iDate},estChInCortexA{iDate}(iRun,:),...
                estChInCortexB{iDate}(iRun,:),gammaRange,lowFreqRange);

            % Get comodulogram
            winSize   = 100e3; 
            stepSize  = 50e3; 
            nHigh     = size(gammaRange,2);
            nLow      = size(lowFreqRange,2);
            nChan     = min(size(lowFreqA,3),size(lowFreqB,3)); 


            % Calculate windows
            winStart = 1:stepSize:(dataLen-winSize);
            winEnd   = winStart+winSize-1; 
            nWin     = numel(winStart); 

            % Preallocating modulation indices
            modA2B = NaN(nWin,nHigh,nLow,nChan); 
            modB2A = NaN(nWin,nHigh,nLow,nChan); 
            modA2A = NaN(nWin,nHigh,nLow,nChan); 
            modB2B = NaN(nWin,nHigh,nLow,nChan);

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

                modA2BT = NaN(nWin,nLow,nChan);
                modB2AT = NaN(nWin,nLow,nChan);
                modA2AT = NaN(nWin,nLow,nChan);
                modB2BT = NaN(nWin,nLow,nChan);

                for iWin = 1:nWin
                    idx = winStart(iWin):winEnd(iWin);
                    for iLow = 1:nLow
                        [modA2BT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)),0);
                        [modB2AT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)),0);
                        [modA2AT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)),0);
                        [modB2BT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)),0);
                    end
                end             
               
                modA2B(:,iHigh,:,:) = modA2BT;
                modB2A(:,iHigh,:,:) = modB2AT;
                modA2A(:,iHigh,:,:) = modA2AT;
                modB2B(:,iHigh,:,:) = modB2BT;

            end
            toc;

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
                'modA2B','modB2A','modA2A','modB2B','-append');
            clear lowFreqAconst lowFreqBconst highFreqAconst highFreqBconst

            modIdxAllA2B{iRun,iDate} = modA2B;
            modIdxAllB2A{iRun,iDate} = modB2A;
            modIdxAllA2A{iRun,iDate} = modA2A;
            modIdxAllB2B{iRun,iDate} = modB2B; 
            toc;        
      
        else
            clear vars;
            vars = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']); 
            modIdxAllA2B{iRun,iDate} = vars.modA2B;
            modIdxAllB2A{iRun,iDate} = vars.modB2A;
            modIdxAllA2A{iRun,iDate} = vars.modA2A;
            modIdxAllB2B{iRun,iDate} = vars.modB2B;

        end
    end
end

%% Getting the shuffled comodulogram
winSize   = 100e3;
stepSize  = 100e3;
nHigh     = size(gammaRange,2);
nLow      = size(lowFreqRange,2);

for iDate = 2:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);

        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);

        if chA(1)== 0 || chB(1)==0; continue; end
%%
        [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
            allProbeData{fileNum,iDate}.probe1Ch,allProbeData{fileNum,iDate}.probe2Ch,badElecA{fileNum,iDate},...
            badElecB{fileNum,iDate},allBadTimes{fileNum,iDate},estChInCortexA{iDate}(iRun,:),...
            estChInCortexB{iDate}(iRun,:),gammaRange,lowFreqRange);
%%
        % Shuffle method 1: Segment the data into different 100 s windows,
        % shuffle the windows, and get the modulation index
        dataLen = size(amplitudeA,2);
        nChan     = min(size(amplitudeA,3),size(amplitudeB,3));

        % Calculate windows
        winStart = 1:stepSize:(dataLen-winSize);
        winEnd   = winStart+winSize-1;
        nWin     = numel(winStart);

        % Generate all shuffle pairs (iWin, jWin) where iWin ~= jWin
        [w1, w2]     = ndgrid(1:nWin, 1:nWin);
        shufflePairs = [w1(:), w2(:)];
        shufflePairs = shufflePairs(shufflePairs(:,1) ~= shufflePairs(:,2), :);
        nShuffle     = size(shufflePairs,1);

        highFreqAconst = parallel.pool.Constant(amplitudeA);
        highFreqBconst = parallel.pool.Constant(amplitudeB);
        lowFreqAconst  = parallel.pool.Constant(phaseA);
        lowFreqBconst  = parallel.pool.Constant(phaseB);

        comb = combvec(1:nHigh, 1:nLow, 1:nShuffle)';
        nComb = size(comb,1);

        tic
        parfor iShuffle = 1: nComb
            iHigh = comb(iShuffle,1);
            iLow  = comb(iShuffle,2);
            iRow  = comb(iShuffle,3);

            idx1 = winStart(shufflePairs(iRow,1)):winEnd(shufflePairs(iRow,1));
            idx2 = winStart(shufflePairs(iRow,2)):winEnd(shufflePairs(iRow,2));

            lowA  = lowFreqAconst.Value;
            lowB  = lowFreqBconst.Value;
            highA = highFreqAconst.Value;
            highB = highFreqBconst.Value;

            modShuffleA2B{iShuffle} = getPhaseAmpCoupling(squeeze(lowB(iLow,idx1,:)),squeeze(highA(iHigh,idx2,:)));
            modShuffleB2A{iShuffle} = getPhaseAmpCoupling(squeeze(lowA(iLow,idx1,:)),squeeze(highB(iHigh,idx2,:)));
            modShuffleA2A{iShuffle} = getPhaseAmpCoupling(squeeze(lowA(iLow,idx1,:)),squeeze(highA(iHigh,idx2,:)));
            modShuffleB2B{iShuffle} = getPhaseAmpCoupling(squeeze(lowB(iLow,idx1,:)),squeeze(highB(iHigh,idx2,:)));
        end
        toc

        for iC = 1:nComb
            iHigh = comb(iC,1);
            iLow  = comb(iC,2);
            iRow  = comb(iC,3);

            modIdxAllA2BShuffleT(iRow,iHigh,iLow,:) = modShuffleA2B{iC};
            modIdxAllB2AShuffleT(iRow,iHigh,iLow,:) = modShuffleB2A{iC};
            modIdxAllA2AShuffleT(iRow,iHigh,iLow,:) = modShuffleA2A{iC};
            modIdxAllB2BShuffleT(iRow,iHigh,iLow,:) = modShuffleB2B{iC};
        end

        modIdxAllA2BShuffle{iRun,iDate} = modIdxAllA2BShuffleT;
        modIdxAllB2AShuffle{iRun,iDate} = modIdxAllB2AShuffleT;
        modIdxAllA2AShuffle{iRun,iDate} = modIdxAllA2AShuffleT;
        modIdxAllB2BShuffle{iRun,iDate} = modIdxAllB2BShuffleT;


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

        combShift  = combvec(1:nHigh, 1:nLow, 1:nShift,1:nWin)';
        nCombShift = size(combShift,1);
        
        tic;
        parfor iShuffle = 1:nCombShift
            iHigh   = combShift(iShuffle,1);
            iLow    = combShift(iShuffle,2);
            iShift  = combShift(iShuffle,3);
            iRow    = combShift(iShuffle,4);

            idx = winStart(iRow):winEnd(iRow);

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

        for iC = 1:nCombShift
            iHigh = combShift(iC,1);
            iLow  = combShift(iC,2);
            iRow  = combShift(iC,3);

            modIdxAllA2CircleT(iRow,iHigh,iLow,:) = modCircleA2B{iC};
            modIdxAllB2ACircleT(iRow,iHigh,iLow,:) = modCircleB2A{iC};
            modIdxAllA2ACircleT(iRow,iHigh,iLow,:) = modCircleA2A{iC};
            modIdxAllB2BCircleT(iRow,iHigh,iLow,:) = modCircleB2B{iC};
        end

        modIdxAllA2Circle{iRun,iDate}  = modIdxAllA2CircleT;
        modIdxAllB2ACircle{iRun,iDate} = modIdxAllB2ACircleT;
        modIdxAllA2ACircle{iRun,iDate} = modIdxAllA2ACircleT;
        modIdxAllB2BCircle{iRun,iDate} = modIdxAllB2BCircleT;
    end
end

%% Getting the recording shuffled comodulogram 



for iDate = 2:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        clear probeA probeB chA chB

        if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')
            varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
            varIdx  = ismember('modShuffleA2B',varInfo);
        else
            varIdx = 1;
        end

        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],'file')|| ~varIdx 
            clc; disp(['Getting control distribution for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

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

            % Get the bandlimited data
            clear highFreqA highFreqB lowFreqA lowFreqB modA2B modB2A

            % High frequencies
            disp('Filtering high frequency signals....')
            for iHigh = 1:size(gammaRange,2)
                fLow = gammaRange(iHigh);
                fHigh = gammaRange(iHigh)+5; % 5 Hz bandwidth
                highFreqA(iHigh,:,:) = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,1e3))';
                highFreqB(iHigh,:,:) = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,1e3))';
            end

            % Low frequencies
            disp('Filtering low frequency signals....')
            for iLow = 1:size(lowFreqRange,2)
                fLow = lowFreqRange(iLow);
                fHigh = lowFreqRange(iLow)+2; % 2 Hz bandwidth
                filtOrder = 3*fs/lowFreqRange(iLow);
                epochLen = round(3*filtOrder);

                if ~mod(dataLen,epochLen)==0
                    if ~isprime(dataLen/1e3)
                        epochLen = ceil(epochLen/1000)*1e3;
                        iCount = 1; 
                        while ~mod(dataLen,epochLen)==0 && iCount~= 5
                            epochLen = epochLen+filtOrder; 
                            iCount = iCount+1; 
                        end
                    else
                        epochLen=0;
                    end
                end

                if iCount == 5
                    epochLen = 0; 
                end 
                try
                    lowFreqA(iLow,:,:) = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,epochLen))';
                    lowFreqB(iLow,:,:) = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,epochLen))';
                catch
                    lowFreqA(iLow,:,:) = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,0))';
                    lowFreqB(iLow,:,:) = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,0))';
                end
            end

            [nLow, ~, nChan]  = size(lowFreqA);
            [nHigh, ~, ~]         = size(highFreqA);
            minShift              = 10e3; % 10 sec shifting
            nSurrogates           = 100;

            % Preallocate
            MISurrA2B = NaN(nSurrogates, nHigh, nLow, nChan, 'single');
            MISurrB2A = NaN(nSurrogates, nHigh, nLow, nChan, 'single');
            MISurrA2A = NaN(nSurrogates, nHigh, nLow, nChan, 'single');
            MISurrB2B = NaN(nSurrogates, nHigh, nLow, nChan, 'single');
          %%  
            tic;
            parfor iRow = 1:nSurrogates
                disp(['Shuffle: ' num2str(iRow)]);
                rng('shuffle');
                shift = randi([minShift dataLen - minShift]); % random shift

                % Shift high-freq analytic signal in time
                highShiftA = circshift(highFreqA,[0,shift, 0]);
                highShiftB = circshift(highFreqB,[0,shift,0]);

                % Compute MI for all freq/channel combinations
                MIA2BT = NaN(nHigh,nLow,nChan,'single');
                MIB2AT = NaN(nHigh,nLow,nChan,'single');
                MIA2AT = NaN(nHigh,nLow,nChan,'single');
                MIB2BT = NaN(nHigh,nLow,nChan,'single');
                
                for iH = 1:nHigh
                    for iL = 1:nLow
                        % Calculate PAC
                        MIA2BT(iH,iL,:) = single(getPhaseAmpCoupling(squeeze(lowFreqB(iL,:,:)),squeeze(highShiftA(iH,:,:)),0));
                        MIB2AT(iH,iL,:) = single(getPhaseAmpCoupling(squeeze(lowFreqA(iL,:,:)),squeeze(highShiftB(iH,:,:)),0));
                        MIA2AT(iH,iL,:) = single(getPhaseAmpCoupling(squeeze(lowFreqA(iL,:,:)),squeeze(highShiftA(iH,:,:)),0));
                        MIB2BT(iH,iL,:) = single(getPhaseAmpCoupling(squeeze(lowFreqB(iL,:,:)),squeeze(highShiftB(iH,:,:)),0));
                    end
                end

                MISurrA2B(iRow,:,:,:) = MIA2BT;
                MISurrB2A(iRow,:,:,:) = MIB2AT;
                MISurrA2A(iRow,:,:,:) = MIA2AT;
                MISurrB2B(iRow,:,:,:) = MIB2BT;
            end
            toc;

            %%
            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
                'MISurrA2B','MISurrB2A','MISurrA2A','MISurrB2B','-append');

            % Get comodulogram
            winSize   = 100e3; 
            stepSize  = 100e3; 
            nHigh     = size(gammaRange,2);
            nLow      = size(lowFreqRange,2);
            nChan     = min(size(lowFreqA,3),size(lowFreqB,3)); 


            % Calculate windows
            winStart = 1:stepSize:(dataLen-winSize);
            winEnd   = winStart+winSize-1; 
            nWin     = numel(winStart); 

            % Avoid overhead
          

            % % Generate all shuffle pairs (iWin, jWin) where iWin ~= jWin
            % [w1, w2]     = ndgrid(1:nWin, 1:nWin);
            % shufflePairs = [w1(:), w2(:)];
            % shufflePairs = shufflePairs(shufflePairs(:,1) ~= shufflePairs(:,2), :);
            % nShuffle     = 15; 
            % rng('shuffle');
            % shuffleRows  = randi([1 size(shufflePairs,1)],1,nShuffle);
            % nShuffle = size(shufflePairs,1);

            modShuffleA2B = NaN(nShuffle,nHigh,nLow,nChan); 
            modShuffleB2A = NaN(nShuffle,nHigh,nLow,nChan);
            modShuffleA2A = NaN(nShuffle,nHigh,nLow,nChan);
            modShuffleB2B = NaN(nShuffle,nHigh,nLow,nChan);

            lowFreqAconst  = parallel.pool.Constant(lowFreqA);
            lowFreqBconst  = parallel.pool.Constant(lowFreqB);
            highShiftA = NaN(nShuffle,nHigh,dataLen,nChan,'single');            
            highShiftB = NaN(nShuffle,nHigh,dataLen,nChan,'single');

            for iSh = 1:nShuffle
                rng('shuffle');
                shift = randi([minShift dataLen - minShift]); % random shift
                highShiftA(iSh,:,:,:) = circshift(highFreqA,[0,shift, 0]);
                highShiftB(iSh,:,:,:) = circshift(highFreqB,[0,shift,0]);
            end
            
            highFreqAconst = parallel.pool.Constant(highShiftA);
            highFreqBconst = parallel.pool.Constant(highShiftB);

            tic;
            comb = combvec(1:nHigh, 1:nLow, 1:15)';
            nComb = size(comb,1);

            clear tempA2B tempB2A tempA2A tempB2B
            % tempA2B = NaN(nWin,nComb,nChan);
            % tempB2A = NaN(nWin,nComb,nChan);
            % tempA2A = NaN(nWin,nComb,nChan);
            % tempB2B = NaN(nWin,nComb,nChan);

%%
tic;
for iC = 1: nComb
    iHigh = comb(iC,1);
    iLow  = comb(iC,2);
    iRow     = comb(iC,3);
    disp(num2str(iC));
    % Local copies
    lowA = lowFreqAconst.Value;
    lowB = lowFreqBconst.Value;
    highA = highFreqAconst.Value;
    highB = highFreqBconst.Value;

    %  for iWin = 1:nWin
    %     idx = winStart(iWin):winEnd(iWin);
    %     for iLow = 1:nLow
    %         [modA2BT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)),0);
    %         [modB2AT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)),0);
    %         [modA2AT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowA(iLow,idx,1:nChan)),squeeze(highA(iHigh,idx,1:nChan)),0);
    %         [modB2BT(iWin,iLow,:),~] = getPhaseAmpCoupling(squeeze(lowB(iLow,idx,1:nChan)),squeeze(highB(iHigh,idx,1:nChan)),0);
    %     end
    % end

    % for iWin = 1: nWin
    rng('shuffle');
    iWin = randi([1 5]);
    idx = winStart(iWin):winEnd(iWin);

    lA = squeeze(lowA(iLow,idx,:));
    lB = squeeze(lowB(iLow,idx,:));
    hA = squeeze(highA(iRow,iHigh,idx,:));
    hB = squeeze(highB(iRow,iHigh,idx,:));

    tempA2B{iC} = getPhaseAmpCoupling(lB, hA, 0);
    tempB2A{iC} = getPhaseAmpCoupling(lA, hB, 0);
    tempA2A{iC} = getPhaseAmpCoupling(lA, hA, 0);
    tempB2B{iC} = getPhaseAmpCoupling(lB, hB, 0);
    % end
end
toc;
%%
            combNew = combvec(1:nHigh, 1:nLow, 1:nShuffle)';
            for iC = 1:nComb
                iHigh = combNew(iC,1);
                iLow  = combNew(iC,2);
                iRow     = combNew(iC,3);

                modShuffleA2B(iRow,iHigh,iLow,:) = tempA2B{iC};% mean(tempA2B,1,'omitnan');
                modShuffleB2A(iRow,iHigh,iLow,:) = tempB2A{iC};%mean(tempB2A,1,'omitnan');
                modShuffleA2A(iRow,iHigh,iLow,:) = tempA2A{iC};%mean(tempA2A,1,'omitnan');
                modShuffleB2B(iRow,iHigh,iLow,:) = tempB2B{iC};%mean(tempB2B,1,'omitnan');
            end

            toc;

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
                'modShuffleA2B','modShuffleB2A','modShuffleA2A','modShuffleB2B','-append');
            clear lowFreqAconst lowFreqBconst highFreqAconst highFreqBconst

            modIdxAllA2BShuffle{iRun,iDate} = modShuffleA2B;
            modIdxAllB2AShuffle{iRun,iDate} = modShuffleB2A;
            modIdxAllA2AShuffle{iRun,iDate} = modShuffleA2A;
            modIdxAllB2BShuffle{iRun,iDate} = modShuffleB2B; 
            toc;

        else
            clear vars;
            vars = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']); 
            modIdxAllA2BShuffle{iRun,iDate} = vars.modShuffleA2B;
            modIdxAllB2AShuffle{iRun,iDate} = vars.modShuffleB2A;
            modIdxAllA2AShuffle{iRun,iDate} = vars.modShuffleA2A;
            modIdxAllB2BShuffle{iRun,iDate} = vars.modShuffleB2B;
        end
    end
end

MIthresh = prctile(x, 99);

%% Plot the comodulogram
avgModA2B = squeeze(mean(modIdxAllA2B{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan')); 
avgModB2A = squeeze(mean(modIdxAllB2A{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
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
avgModA2A = squeeze(mean(modIdxAllA2A{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllA2B{5,2},1,'omitnan')); 
avgModB2B = squeeze(mean(modIdxAllB2B{iRun,iDate} ,1,'omitnan'));%squeeze(mean(modIdxAllB2A{5,2},1,'omitnan'));
figure; contourf(lowFreqRange,gammaRange,median(avgModA2A,3,'omitnan'),'lines','none'); clim([0 1e-4]); colorbar; colormap jet;
figure; contourf(lowFreqRange,gammaRange,median(avgModB2B,3,'omitnan'),'lines','none');clim([0 1e-4]) ;colorbar; colormap jet;
 


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
    varIdx  = ismember('amplitudeA',varInfo);
else
    varIdx = 1;
end

if ~varIdx
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

    % High frequencies
    disp('Getting high frequency amplitudes ....');
    for iHigh = 1:size(gammaRange,2)
        fLow = gammaRange(iHigh);
        fHigh = gammaRange(iHigh)+5; % 5 Hz bandwidth
        highFreqA(iHigh,:,:) = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,1e3))';
        highFreqB(iHigh,:,:) = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,1e3))';
    end
    amplitudeA = abs(hilbert(highFreqA));
    amplitudeB = abs(hilbert(highFreqB));

    % Low frequencies
    disp('Filtering low frequency signals....');
    for iLow = 1:size(lowFreqRange,2)
        fLow = lowFreqRange(iLow);
        fHigh = lowFreqRange(iLow)+2; % 2 Hz bandwidth
        filtOrder = 3*fs/lowFreqRange(iLow);
        epochLen = round(3*filtOrder);

        if ~mod(dataLen,epochLen)==0
            if ~isprime(dataLen/1e3)
                epochLen = ceil(epochLen/1000)*1e3;
                iCount = 1;
                while ~mod(dataLen,epochLen)==0 && iCount~= 5
                    epochLen = epochLen+filtOrder;
                    iCount = iCount+1;
                end
            else
                epochLen=0;
            end
        end

        if iCount == 5
            epochLen = 0;
        end

        try
            lowFreqA(iLow,:,:)  = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,epochLen))';
            lowFreqB(iLow,:,:)  = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,epochLen))';

        catch
            lowFreqA(iLow,:,:)  = single(eegfilt(probeA(1:dataLen,chA(1):chA(2))',fs,fLow,fHigh,0))';
            lowFreqB(iLow,:,:)  = single(eegfilt(probeB(1:dataLen,chB(1):chB(2))',fs,fLow,fHigh,0))';
        end
    end

    phaseA = rad2deg(angle(hilbert(lowFreqA))); % Get phase of low frequency
    phaseB = rad2deg(angle(hilbert(lowFreqB))); % Get phase of low frequency


    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat'],...
        'amplitudeA','amplitudeB','phaseA','phaseB','-append');
else
    vars = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_' num2str(fileNum) '.mat']);
    amplitudeA = vars.amplitudeA;
    amplitudeB = vars.amplitudeB;
    phaseA     = vars.phaseA;
    phaseB     = vars.phaseB;
end
end