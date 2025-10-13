%% ephysDualProbeRS_v6
% This function performs analysis on LFP recorded simultaneously from two linear electrode arrays.
% June 25,2025 - Keerthana Manikandan
% This code performs the following for ONE MONKEY (for combined analyis,
% check this script- combinedAnalysisEphysDualProbeRS.m)
% Check ephysDualProbeRS_v5 for previous iterations

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
saveFigureFlag = 1;


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

%% Store/Retrieve LFP from the two probes...
disp(['Storing/Retrieving LFP data for ' monkeyName]);

for iDate = 1: size(allDates,1)
    clear expDate
    expDate = allDates(iDate,:);
    saveFolder = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];

    for iRun = 1:length(datFileNumAll{iDate,1})
        clear fileNum zL pL kL
        fileNum = datFileNumAll{iDate,1}(iRun);

        % Get the name of stored file
        if strcmp(expDate,'11_01_2021') && (fileNum == 1 || fileNum == 2) % Only for Charlie Sheen
            datFileName = 'datafile_000';
        else
            datFileName = datFileNameAll{iDate,1};
        end

        if (fileNum>=10)
            datFileName = datFileName(1:end-1);
        end

        % Check if the LFP is already stored or not
        if  ~exist([saveFolder '\' datFileName num2str(datFileNumAll{iDate,1}(iRun)) '_lfp.mat'],'file')

            [probe1Ch,probe2Ch,raw1Ch,raw2Ch,eegCh,scalpEEGCh] = ...
                saveLFPDualProbe(monkeyName,expDate,fileNum,datFileName,serverPath,fs);

            % Save the LFP
            if ~exist('saveFolder','dir'); [~,~] = mkdir(saveFolder); end
            disp('Storing data... ' );

            if exist('eegCh','var') && exist('scalpEEGCh','var')
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch','raw1Ch','raw2Ch','eegCh','scalpEEGCh');

            elseif exist('eegCh','var') && ~exist('scalpEEGCh','var')
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','raw1Ch','raw2Ch','probe2Ch','eegCh');

            elseif  ~(exist('eegCh','var') || exist('scalpEEGCh','var'))
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch','raw1Ch','raw2Ch');
            end

            allProbeData{fileNum,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);

        else
            % Retrieve LFP Data
            disp(['Retrieving data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            allProbeData{fileNum,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
        end
    end
end

clear probe1Ch probe2Ch eegCh
clc; disp('Data Stored/Retrieved');

% %% Convert doubles to singles
% for iDate = 1:size(allDates,1)
%     clear expDate datFileNum
%     expDate    = allDates(iDate,:);
%     datFileNum = datFileNumAll{iDate,1};
%     for iRun = 1: length(datFileNum)
%         saveFolder = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
%         fileNum    = datFileNumAll{iDate,1}(iRun);
%
%         if strcmp(expDate,'11_01_2021') && (fileNum == 1 || fileNum == 2) % Only for Charlie Sheen
%             datFileName = 'datafile_000';
%         else
%             datFileName = datFileNameAll{iDate,1};
%         end
%
%         if (fileNum>=10)
%             datFileName = datFileName(1:end-1);
%         end
%
%         disp([ expDate ' File: ' num2str(fileNum)]);
%
%         probe1Ch   = single(allProbeData{fileNum,iDate}.probe1Ch);
%         probe2Ch   = single(allProbeData{fileNum,iDate}.probe2Ch);
%         raw1Ch     = single(allProbeData{fileNum,iDate}.raw1Ch);
%         raw2Ch     = single(allProbeData{fileNum,iDate}.raw2Ch);
%         eegCh      = single(allProbeData{fileNum,iDate}.eegCh);
%         scalpEEGCh = single(allProbeData{fileNum,iDate}.scalpEEGCh);
%
%         if exist('eegCh','var') && exist('scalpEEGCh','var')
%             save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch','raw1Ch','raw2Ch','eegCh','scalpEEGCh');
%
%         elseif exist('eegCh','var') && ~exist('scalpEEGCh','var')
%             save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','raw1Ch','raw2Ch','probe2Ch','eegCh');
%
%         elseif  ~(exist('eegCh','var') || exist('scalpEEGCh','var'))
%             save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch','raw1Ch','raw2Ch');
%         end
%
%         allProbeData{fileNum,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
%
%
%     end
% end

%% Determine bad time segments and channels from the two probes
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\badChTimeVars.mat'],'file') 
    tic;
    [allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,allProbeData,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);
    toc;
    disp(['Identified bad channels and time segments for ' monkeyName]);
    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\badChTimeVars.mat'],'allBadTimes','badElecA','badElecB');
else
    load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\badChTimeVars.mat']);
end

% Get the transition channels for the recordings...
% Eliminate bad channels and obtain the transition channels by
% computing slope of the marginals (obtained by averaging the intra-probe
% correlograms) - Mean is used here for sensitivity to outliers which is
% needed for picking the transition channel
clear estChInCortexA estChInCortexB
for iDate = 1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    clc; disp(expDate);

    for iRun = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iRun == 4;  estChInCortexA{iDate}(iRun,:) = 0;  estChInCortexB{iDate}(iRun,:) = 0; continue; end
        if strcmp(expDate,'02_07_2023') && iRun == 2;  estChInCortexA{iDate}(iRun,:) = 0;  estChInCortexB{iDate}(iRun,:) = 0; continue; end
        if strcmp(expDate,'04_11_2023') && iRun == 10; estChInCortexA{iDate}(iRun,:) = 0;  estChInCortexB{iDate}(iRun,:) = 0; continue; end
        if strcmp(expDate,'11_01_2021') && iRun == 9;  estChInCortexA{iDate}(iRun,:) = 0;  estChInCortexB{iDate}(iRun,:) = 0; continue; end
        if strcmp(expDate,'01_11_2022') && iRun == 6;  estChInCortexA{iDate}(iRun,:) = 0;  estChInCortexB{iDate}(iRun,:) = 0; continue; end 
        fileNum = datFileNum(iRun);

        if isempty(allProbeData{fileNum,iDate}.probe1Ch) || isempty(allProbeData{fileNum,iDate}.probe2Ch); continue; end % Check if probes do not have any LFP

        % Get the ephys data
        clear probeA probeB marginalVal fx
        probeA = allProbeData{fileNum,iDate}.probe1Ch;
        probeB = allProbeData{fileNum,iDate}.probe2Ch;

        % Remove bad channels
        probeA(:,badElecA{fileNum,iDate}) = [];
        probeB(:,badElecB{fileNum,iDate}) = [];

        % Remove bad times
        probeA(allBadTimes{fileNum,iDate},:) = [];
        probeB(allBadTimes{fileNum,iDate},:) = [];

        % Calculate the envelope/power for the LFPs
        probeA = envelope(probeA,5);
        probeB = envelope(probeB,5);

        % Get the transition channel as determined during the experiment
        transitionChA = chInCortexProbeA{iDate}(iRun);
        transitionChB = chInCortexProbeB{iDate}(iRun);

        for iProbe = 1:2
            clear probeTemp transitionCh estChInCortex
            switch iProbe
                case 1 % Probe A
                    probeTemp = probeA;
                    transitionCh = transitionChA;

                case 2 % Probe B
                    probeTemp = probeB;
                    transitionCh = transitionChB;
            end

            % Determine the channels that are within cortex
            if size(probeTemp,2)~=1

                % Get mean within probe marginals from wideband range for the
                % electrode
                marginalVal = mean(imgaussfilt(corr(probeTemp,'Rows','complete'),1),2,'omitnan');

                % Get the slope of the marginals for the electrode
                fx = abs(movmean(gradient(marginalVal),2,'omitnan'));

                % Find the channel with maximum slope
                clear tempCh
                tempCh = find(fx == max(fx));

                if tempCh>=20 % Set transition to 1 if max slope is after channel 20
                    tempCh = 1;
                end

                if abs(tempCh - transitionCh)>=10 % Set transition to the value determined from notes if difference between estimated and observed value exceeds 10
                    estChInCortex(1) = transitionCh;

                elseif abs(tempCh - transitionCh)<= 3 % Set transition to estimated value if difference between estimated and observed value is less than or equal to 3
                    estChInCortex(1) = tempCh;
                else
                    estChInCortex(1) = floor((tempCh + transitionCh)/2); % Set transition to mean between estimated and observed value otherwise
                end

                if estChInCortex+ 20 > size(probeTemp,2) % Limit to the first 20 channels only
                    estChInCortex(2)= size(probeTemp,2);
                else
                    estChInCortex(2)= estChInCortex+ 20;
                end 

            else
                estChInCortex = [1 1];
            end
            % Check the difference between estimated and observed
            % transition
            diffCalcObs(iDate,iRun) = abs(transitionCh-estChInCortex(1));
            switch iProbe
                case 1
                    estChInCortexA{iDate}(iRun,:) = estChInCortex;
                    diffCalcObsA = diffCalcObs;
                case 2
                    estChInCortexB{iDate}(iRun,:) = estChInCortex;
                    diffCalcObsB = diffCalcObs;
            end
        end
    end
end
disp(['Determined transition channels for ' monkeyName]);

%% Separate LFP into different frequencies and correlate between the electrodes (within frequency)
clear meanPSDEEG medPSDEEG  specAMeanAll specBMeanAll specAMedAll specBMedAll meanCorrInfraSlow medCorrInfraSlow
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

bandLabels = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking'};
timeLabels = {'Time series','Power','Infraslow'};

% Check if certain variables are stored....
if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') 
    varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    varIdx  = ismember('infraIntraAAllR',varInfo);
else 
    varIdx = 1;
end

% Calculate pairwise correlations between and within probes
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') || ~varIdx
    tic; % re-run for Whiskey... 
    for iDate = 1: size(allDates,1)
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

            rawA   = allProbeData{fileNum,iDate}.raw1Ch;
            rawB   = allProbeData{fileNum,iDate}.raw2Ch;

            % Remove bad channels
            probeA(:,badElecA{fileNum,iDate}) = [];
            probeB(:,badElecB{fileNum,iDate}) = [];

            rawA(:,badElecA{fileNum,iDate}) = [];
            rawB(:,badElecB{fileNum,iDate}) = [];

            % Remove bad times
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];

            rawA(allBadTimes{fileNum,iDate},:) = [];
            rawB(allBadTimes{fileNum,iDate},:) = [];

            chA = estChInCortexA{iDate}(iRun,:);
            chB = estChInCortexB{iDate}(iRun,:);

            if chA(1) == 0 || chB(1) == 0 || isempty(probeA) || isempty(probeB)                
                meanPairCorr(iDate,iRun,1:5)   = NaN; maxPairCorr(iDate,iRun,1:5)    = NaN;
                medPairCorr(iDate,iRun,1:5)    = NaN; meanIntraCorrA(iDate,iRun,1:5) = NaN;
                meanIntraCorrB(iDate,iRun,1:5) = NaN; medIntraCorrA(iDate,iRun,1:5)  = NaN;
                medIntraCorrB(iDate,iRun,1:5)  = NaN;

                medPairCorrSuper(iDate,iRun,1:5)      = NaN;
                medPairCorrMid(iDate,iRun,1:5)        = NaN;
                medPairCorrDeep(iDate,iRun,1:5)       = NaN;

                envelopePairCorrSuper(iDate,iRun,1:5) = NaN;
                envelopePairCorrMid(iDate,iRun,1:5)   = NaN;
                envelopePairCorrDeep(iDate,iRun,1:5)  = NaN;

                infraPairCorrSuper(iDate,iRun,1:5)    = NaN;
                infraPairCorrMid(iDate,iRun,1:5)      = NaN;
                infraPairCorrDeep(iDate,iRun,1:5)     = NaN;

                continue;
            end

            % Get mean, median and maximum pairwise correlations for different
            % frequency bands...
            clear timeSeriesCorr sizeSpecA sizeSpecB specA_R specB_R freqSpecCorr
            for iBand = 1:5
                clear xA xB xARef xBRef specA specB specEEG timeValsSpec freqValsSpec eegGood envelopeABandLimited envelopeBBandLimited infraSlowA infraSlowB

                if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 ))||...
                        (strcmp(expDate,'02_07_2023') && (fileNum == 2 ))... %||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))...
                        ||(strcmp(expDate,'09_19_2022') && (fileNum == 4)) || (isempty(probeA) && isempty(probeB))...
                        || (strcmp(expDate,'04_11_2023') && (fileNum == 10))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||

                    meanPairCorr(iDate,iRun,iBand)   = NaN; maxPairCorr(iDate,iRun,iBand)    = NaN; medPairCorr(iDate,iRun,iBand)  = NaN;
                    meanIntraCorrA(iDate,iRun,iBand) = NaN; meanIntraCorrB(iDate,iRun,iBand) = NaN;
                    medIntraCorrA(iDate,iRun,iBand)  = NaN; medIntraCorrB(iDate,iRun,iBand)  = NaN;

                    medPairCorrSuper(iDate,iRun,iBand)  = NaN; medPairCorrMid(iDate,iRun,iBand)  = NaN; medPairCorrDeep(iDate,iRun,iBand)  = NaN;
                    envelopePairCorrSuper(iDate,iRun,iBand)  = NaN; envelopePairCorrMid(iDate,iRun,iBand)  = NaN; envelopePairCorrDeep(iDate,iRun,iBand)  = NaN;
                    infraPairCorrSuper(iDate,iRun,iBand)  = NaN; infraPairCorrMid(iDate,iRun,iBand)  = NaN; infraPairCorrDeep(iDate,iRun,iBand)  = NaN;

                    continue;
                end

                if chA(1) == 0 || chB(1)== 0 % Check if any of the recordings are empty
                    meanPairCorr(iDate,iRun,iBand)   = NaN; maxPairCorr(iDate,iRun,iBand)    = NaN;
                    medPairCorr(iDate,iRun,iBand)    = NaN; meanIntraCorrA(iDate,iRun,iBand) = NaN;
                    meanIntraCorrB(iDate,iRun,iBand) = NaN; medIntraCorrA(iDate,iRun,iBand)  = NaN;
                    medIntraCorrB(iDate,iRun,iBand)  = NaN;
                    continue;
                end

                clear xA yA
                switch iBand % Get the correlation vs connectivity vs distance plots for different bands...
                    case 1 % Theta band
                        xA = filtfilt(bT,aT,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bT,aT,double(probeB(:,chB(1):chB(2))));

                    case 2 % Alpha band
                        xA = filtfilt(bA,aA,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bA,aA,double(probeB(:,chB(1):chB(2))));

                    case 3 % Beta band
                        xA = filtfilt(bB,aB,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bB,aB,double(probeB(:,chB(1):chB(2))));

                    case 4 % Gamma band
                        xA = filtfilt(bG,aG,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bG,aG,double(probeB(:,chB(1):chB(2))));

                    case 5 % Spiking
                        xA = rawA(:,chA(1):chA(2));
                        xB = rawB(:,chB(1):chB(2));
                end

                % Get the intra probe correlations for channels inside the
                % cortex...
                medIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(xA,'rows','complete'))),'all','omitnan');
                medIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(xB,'rows','complete'))),'all','omitnan');
               
                intraCorrA{iDate,iRun,iBand} = single(tril(corr(xA,'rows','complete')));
                intraCorrB{iDate,iRun,iBand} = single(tril(corr(xB,'rows','complete')));

                % Get pairwise correlations between the two probes...
                maxPairCorr(iDate,iRun,iBand)  = max(corr(xA,xB),[],'all','omitnan');
                meanPairCorr(iDate,iRun,iBand) = mean(corr(xA,xB),'all','omitnan');
                medPairCorr(iDate,iRun,iBand)  = median(corr(xA,xB),'all','omitnan');%

                % Get instantaneous power and correlate the powers
                envelopeABandLimited = envelope(abs(xA),5);
                envelopeBBandLimited = envelope(abs(xB),5);

                % Within probe correlations - Envelope
                envelopeIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(envelopeABandLimited,envelopeABandLimited))),'all','omitnan');
                envelopeIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(envelopeBBandLimited,envelopeBBandLimited))),'all','omitnan');

                envelopeIntraAAll{iDate,iRun,iBand} =  single(tril(corr(envelopeABandLimited,'rows','complete')));
                envelopeIntraBAll{iDate,iRun,iBand} =  single(tril(corr(envelopeBBandLimited,'rows','complete')));

                % Correlate instantaneous band power
                medCorrEnvelope(iDate,iRun,iBand)  =  median(corr(envelopeABandLimited,envelopeBBandLimited),'all','omitnan');

                enSizeA = size(envelopeABandLimited);
                enSizeB = size(envelopeBBandLimited);

                % Filter envelope from 0.01-0.1 Hz
                envelopeABandLimited = [envelopeABandLimited; envelopeABandLimited ;envelopeABandLimited];
                infraSlowA = filtfilt(sos,g,double(envelopeABandLimited));
                infraSlowA = single(infraSlowA(enSizeA(1)+1:(end-enSizeA(1)),:));

                envelopeBBandLimited = [envelopeBBandLimited; envelopeBBandLimited ;envelopeBBandLimited];
                infraSlowB = filtfilt(sos,g,double(envelopeBBandLimited));
                infraSlowB = single(infraSlowB(enSizeB(1)+1:(end-enSizeB(1)),:));

                % Within probe correlations - Infraslow
                infraIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(infraSlowA,infraSlowA))),'all','omitnan');
                infraIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(infraSlowB,infraSlowB))),'all','omitnan');

                infraIntraAAll{iDate,iRun,iBand} =  single(tril(corr(infraSlowA,'rows','complete')));
                infraIntraBAll{iDate,iRun,iBand} =  single(tril(corr(infraSlowB,'rows','complete')));

                % Correlate infraslow flucutuations in instantaneous band power
                medCorrInfraSlow(iDate,iRun,iBand)  =  median(corr(infraSlowA,infraSlowB),'all','omitnan');


                % Split the data into superficial, middle or deep for both channels
                chSplit = 6;
                % Electrode A
                if chA(2)-chA(1)== 0 % Single channel
                    timeASuper = xA;
                    timeAMid   = xA;
                    timeADeep  = xA;

                    envelopeASuper = envelopeABandLimited;
                    envelopeAMid   = envelopeABandLimited;
                    envelopeADeep  = envelopeABandLimited;

                    infraASuper = infraSlowA;
                    infraAMid   = infraSlowA;
                    infraADeep  = infraSlowA;
                else
                    timeASuper = xA(:,1:chSplit);
                    timeAMid   = xA(:,(chSplit+1: chSplit*2));
                    timeADeep  = xA(:,2*chSplit+1:end);

                    envelopeASuper = envelopeABandLimited(:,1:chSplit);
                    envelopeAMid   = envelopeABandLimited(:,(chSplit+1: chSplit*2));
                    envelopeADeep  = envelopeABandLimited(:,2*chSplit+1:end);

                    infraASuper = infraSlowA(:,1:chSplit);
                    infraAMid   = infraSlowA(:,(chSplit+1:chSplit*2));
                    infraADeep  = infraSlowA(:,2*chSplit+1:end);

                end

                % Electrode B
                if chB(2)-chB(1)== 0 % Single channel
                    timeBSuper = xB;
                    timeBMid   = xB;
                    timeBDeep  = xB;

                    envelopeBSuper = envelopeBBandLimited;
                    envelopeBMid   = envelopeBBandLimited;
                    envelopeBDeep  = envelopeBBandLimited;

                    infraBSuper = infraSlowB;
                    infraBMid   = infraSlowB;
                    infraBDeep  = infraSlowB;
                else
                    timeBSuper = xB(:,1:chSplit);%
                    timeBMid   = xB(:,(chSplit+1: chSplit*2));%
                    timeBDeep  = xB(:,2*chSplit+1:end);%

                    envelopeBSuper = envelopeBBandLimited(:,1:chSplit);%
                    envelopeBMid   = envelopeBBandLimited(:,(chSplit+1: chSplit*2));%
                    envelopeBDeep  = envelopeBBandLimited(:,2*chSplit+1:end);%

                    infraBSuper = infraSlowB(:,1:chSplit); %
                    infraBMid   = infraSlowB(:,(chSplit+1: chSplit*2));%
                    infraBDeep  = infraSlowB(:,2*chSplit+1:end);%
                end

                % Correlate within compartments
                medPairCorrSuper(iDate,iRun,iBand) = median(corr(timeASuper,timeBSuper,'rows','complete'),'all','omitnan');
                medPairCorrMid(iDate,iRun,iBand)   = median(corr(timeAMid,timeBMid,'rows','complete'),'all','omitnan');
                medPairCorrDeep(iDate,iRun,iBand)  = median(corr(timeADeep,timeBDeep,'rows','complete'),'all','omitnan');

                envelopePairCorrSuper(iDate,iRun,iBand) = median(corr(envelopeASuper,envelopeBSuper,'rows','complete'),'all','omitnan');
                envelopePairCorrMid(iDate,iRun,iBand)   = median(corr(envelopeAMid,envelopeBMid,'rows','complete'),'all','omitnan');
                envelopePairCorrDeep(iDate,iRun,iBand)  = median(corr(envelopeADeep,envelopeBDeep,'rows','complete'),'all','omitnan');

                infraPairCorrSuper(iDate,iRun,iBand) = median(corr(infraASuper,infraBSuper,'rows','complete'),'all','omitnan');
                infraPairCorrMid(iDate,iRun,iBand)   = median(corr(infraAMid,infraBMid,'rows','complete'),'all','omitnan');
                infraPairCorrDeep(iDate,iRun,iBand)  = median(corr(infraADeep,infraBDeep,'rows','complete'),'all','omitnan');

                % Correlate between compartments between probes
                medPairASuperBMid(iDate,iRun,iBand)  = median(corr(timeASuper,timeBMid,'rows','complete'),'all','omitnan');
                medPairASuperBDeep(iDate,iRun,iBand) = median(corr(timeASuper,timeBDeep,'rows','complete'),'all','omitnan');
                medPairAMidBSuper(iDate,iRun,iBand)  = median(corr(timeAMid,timeBSuper,'rows','complete'),'all','omitnan');
                medPairAMidBDeep(iDate,iRun,iBand)   = median(corr(timeAMid,timeBDeep,'rows','complete'),'all','omitnan');
                medPairADeepBSuper(iDate,iRun,iBand) = median(corr(timeADeep,timeBSuper,'rows','complete'),'all','omitnan');
                medPairADeepBMid(iDate,iRun,iBand)   = median(corr(timeADeep,timeBMid,'rows','complete'),'all','omitnan');

                envelopeASuperBMid(iDate,iRun,iBand)  = median(corr(envelopeASuper,envelopeBMid,'rows','complete'),'all','omitnan');
                envelopeASuperBDeep(iDate,iRun,iBand) = median(corr(envelopeASuper,envelopeBDeep,'rows','complete'),'all','omitnan');
                envelopeAMidBSuper(iDate,iRun,iBand)  = median(corr(envelopeAMid,envelopeBSuper,'rows','complete'),'all','omitnan');
                envelopeAMidBDeep(iDate,iRun,iBand)   = median(corr(envelopeAMid,envelopeBDeep,'rows','complete'),'all','omitnan');
                envelopeADeepBSuper(iDate,iRun,iBand) = median(corr(envelopeADeep,envelopeBSuper,'rows','complete'),'all','omitnan');
                envelopeADeepBMid(iDate,iRun,iBand)   = median(corr(envelopeADeep,envelopeBMid,'rows','complete'),'all','omitnan');

                infraASuperBMid(iDate,iRun,iBand)  = median(corr(infraASuper,infraBMid,'rows','complete'),'all','omitnan');
                infraASuperBDeep(iDate,iRun,iBand) = median(corr(infraASuper,infraBDeep,'rows','complete'),'all','omitnan');
                infraAMidBSuper(iDate,iRun,iBand)  = median(corr(infraAMid,infraBSuper,'rows','complete'),'all','omitnan');
                infraAMidBDeep(iDate,iRun,iBand)   = median(corr(infraAMid,infraBDeep,'rows','complete'),'all','omitnan');
                infraADeepBSuper(iDate,iRun,iBand) = median(corr(infraADeep,infraBSuper,'rows','complete'),'all','omitnan');
                infraADeepBMid(iDate,iRun,iBand)   = median(corr(infraADeep,infraBMid,'rows','complete'),'all','omitnan');

            end
        end
    end
    

    % Reshape variables for storage
    disp('Reshaping all variables to store as mat files...')
    matSize   = size(connValsAll);
    connValsR = reshape(connValsAll',[matSize(1)*matSize(2) 1]);
    distValsR = reshape(distSitesAll',[matSize(1)*matSize(2) 1]);

    intraCorrAR    = reshape(intraCorrA,[matSize(1)*matSize(2) size(intraCorrA,3)]);
    intraCorrBR    = reshape(intraCorrB,[matSize(1)*matSize(2) size(intraCorrB,3)]);

    envelopeIntraAAllR = reshape(envelopeIntraAAll,[matSize(1)*matSize(2) size(intraCorrA,3)]);
    envelopeIntraBAllR = reshape(envelopeIntraBAll,[matSize(1)*matSize(2) size(intraCorrA,3)]);

    infraIntraAAllR  = reshape(infraIntraAAll,[matSize(1)*matSize(2) size(intraCorrA,3)]); 
    infraIntraBAllR  = reshape(infraIntraAAll,[matSize(1)*matSize(2) size(intraCorrA,3)]);

    medPairCorrR      = reshape(medPairCorr,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medCorrEnvelopeR  = reshape(medCorrEnvelope,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    medCorrInfraSlowR = reshape(medCorrInfraSlow,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    medIntraCorrAR    = reshape(medIntraCorrA,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    medIntraCorrBR    = reshape(medIntraCorrB,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    envelopeIntraCorrAR = reshape(envelopeIntraCorrA,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    envelopeIntraCorrBR = reshape(envelopeIntraCorrB,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    infraIntraCorrAR = reshape(infraIntraCorrA,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    infraIntraCorrBR = reshape(infraIntraCorrB,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    medPairCorrSuperR = reshape(medPairCorrSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairCorrMidR   = reshape(medPairCorrMid,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    medPairCorrDeepR  = reshape(medPairCorrDeep,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    envelopePairCorrSuperR = reshape(envelopePairCorrSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopePairCorrMidR   = reshape(envelopePairCorrMid,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    envelopePairCorrDeepR  = reshape(envelopePairCorrDeep,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    infraPairCorrSuperR = reshape(infraPairCorrSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraPairCorrMidR   = reshape(infraPairCorrMid,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    infraPairCorrDeepR  = reshape(infraPairCorrDeep,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    medPairASuperBMidR  = reshape(medPairASuperBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairASuperBDeepR = reshape(medPairASuperBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairAMidBSuperR  = reshape(medPairAMidBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairAMidBDeepR   = reshape(medPairAMidBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairADeepBSuperR = reshape(medPairADeepBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairADeepBMidR   = reshape(medPairADeepBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);

    envelopeASuperBMidR  = reshape(envelopeASuperBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeASuperBDeepR = reshape(envelopeASuperBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeAMidBSuperR  = reshape(envelopeAMidBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeAMidBDeepR   = reshape(envelopeAMidBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeADeepBSuperR = reshape(envelopeADeepBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeADeepBMidR   = reshape(envelopeADeepBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);

    infraASuperBMidR  = reshape(infraASuperBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraASuperBDeepR = reshape(infraASuperBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraAMidBSuperR  = reshape(infraAMidBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraAMidBDeepR   = reshape(infraAMidBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraADeepBSuperR = reshape(infraADeepBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraADeepBMidR   = reshape(infraADeepBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);

    nanVals       = (isnan(connValsR) | isnan(medPairCorrR(:,1)));
    lowCorrVals   = (medIntraCorrAR(:,4)<=0.2 | medIntraCorrBR(:,4)<=0.2);
    removeDataIdx = nanVals | lowCorrVals;

    connValsR(removeDataIdx)             = [];
    distValsR(removeDataIdx)             = [];
    medPairCorrR(removeDataIdx,:)        = [];
    medCorrEnvelopeR(removeDataIdx,:)    = [];
    medCorrInfraSlowR(removeDataIdx,:)   = [];
    medIntraCorrAR(removeDataIdx,:)      = [];
    medIntraCorrBR(removeDataIdx,:)      = [];
    envelopeIntraCorrAR(removeDataIdx,:) = [];
    envelopeIntraCorrBR(removeDataIdx,:) = [];
    infraIntraCorrAR(removeDataIdx,:)    = [];
    infraIntraCorrBR(removeDataIdx,:)    = [];

    medPairCorrSuperR(removeDataIdx,:)      = [];
    medPairCorrMidR(removeDataIdx,:)        = [];
    medPairCorrDeepR(removeDataIdx,:)       = [];
    envelopePairCorrSuperR(removeDataIdx,:) = [];
    envelopePairCorrMidR(removeDataIdx,:)   = [];
    envelopePairCorrDeepR(removeDataIdx,:)  = [];
    infraPairCorrSuperR(removeDataIdx,:)    = [];
    infraPairCorrMidR(removeDataIdx,:)      = [];
    infraPairCorrDeepR(removeDataIdx,:)     = [];

    medPairASuperBMidR(removeDataIdx,:)  = [];
    medPairASuperBDeepR(removeDataIdx,:) = [];
    medPairAMidBSuperR(removeDataIdx,:)  = [];
    medPairAMidBDeepR(removeDataIdx,:)   = [];
    medPairADeepBSuperR(removeDataIdx,:) = [];
    medPairADeepBMidR(removeDataIdx,:)   = [];

    envelopeASuperBMidR(removeDataIdx,:)  = [];
    envelopeASuperBDeepR(removeDataIdx,:) = [];
    envelopeAMidBSuperR(removeDataIdx,:)  = [];
    envelopeAMidBDeepR(removeDataIdx,:)   = [];
    envelopeADeepBSuperR(removeDataIdx,:) = [];
    envelopeADeepBMidR(removeDataIdx,:)   = [];

    infraASuperBMidR(removeDataIdx,:)  = [];
    infraASuperBDeepR(removeDataIdx,:) = [];
    infraAMidBSuperR(removeDataIdx,:)  = [];
    infraAMidBDeepR(removeDataIdx,:)   = [];
    infraADeepBSuperR(removeDataIdx,:) = [];
    infraADeepBMidR(removeDataIdx,:)   = [];

    intraCorrAR(removeDataIdx,:) = [];    
    intraCorrBR(removeDataIdx,:) = [];  

    envelopeIntraAAllR(removeDataIdx,:) = []; 
    envelopeIntraBAllR(removeDataIdx,:) = []; 

    infraIntraAAllR(removeDataIdx,:) = []; 
    infraIntraBAllR(removeDataIdx,:) = [];


    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'connValsR','distValsR',...
        'medPairCorrR','medCorrEnvelopeR','medCorrInfraSlowR','medIntraCorrAR','medIntraCorrBR','envelopeIntraCorrAR',...
        'envelopeIntraCorrBR','infraIntraCorrAR','infraIntraCorrBR','medPairCorrSuperR','medPairCorrMidR',...
        'medPairCorrDeepR','envelopePairCorrSuperR','envelopePairCorrMidR','envelopePairCorrDeepR',...
        'infraPairCorrSuperR','infraPairCorrMidR','infraPairCorrDeepR','medPairASuperBMidR','medPairASuperBDeepR',...
        'medPairAMidBSuperR','medPairAMidBDeepR','medPairADeepBSuperR','medPairADeepBMidR','envelopeASuperBMidR',...
        'envelopeASuperBDeepR','envelopeAMidBSuperR','envelopeAMidBDeepR','envelopeADeepBSuperR','envelopeADeepBMidR',...
        'infraASuperBMidR','infraASuperBDeepR','infraAMidBSuperR','infraAMidBDeepR','infraADeepBSuperR','infraADeepBMidR',...
        'removeDataIdx','intraCorrAR','intraCorrBR','envelopeIntraAAllR','envelopeIntraBAllR','infraIntraAAllR','infraIntraBAllR','-append');

toc;
else
    disp('Loading saved variables...')
    clear allVars
    allVars = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    fieldNames = fieldnames(allVars); 

    for iL = 1:length(fieldNames)
        nameVal = fieldNames{iL};
        assignin('base',nameVal,allVars.(nameVal));
    end
    clear allVars;
end

% Check if bad runs are completely removed from the processed data
goodRunsR                   = reshape(goodRuns',[numel(goodRuns) 1]);

% Remove recordings with a single channel
singleChRow = cellfun(@(x) isscalar(x),intraCorrBR(:,1));

goodRunsR(isnan(goodRunsR)) = 0;
goodRunsR(removeDataIdx,:)  = [];


connValsR(~goodRunsR)             = [];
distValsR(~goodRunsR)             = [];
medPairCorrR(~goodRunsR,:)        = [];
medCorrEnvelopeR(~goodRunsR,:)    = [];
medCorrInfraSlowR(~goodRunsR,:)   = [];
medIntraCorrAR(~goodRunsR,:)      = [];
medIntraCorrBR(~goodRunsR,:)      = [];
envelopeIntraCorrAR(~goodRunsR,:) = [];
envelopeIntraCorrBR(~goodRunsR,:) = [];
infraIntraCorrAR(~goodRunsR,:)    = [];
infraIntraCorrBR(~goodRunsR,:)    = [];

medPairCorrSuperR(~goodRunsR|singleChRow,:)      = [];
medPairCorrMidR(~goodRunsR|singleChRow,:)        = [];
medPairCorrDeepR(~goodRunsR|singleChRow,:)       = [];
envelopePairCorrSuperR(~goodRunsR|singleChRow,:) = [];
envelopePairCorrMidR(~goodRunsR|singleChRow,:)   = [];
envelopePairCorrDeepR(~goodRunsR|singleChRow,:)  = [];
infraPairCorrSuperR(~goodRunsR|singleChRow,:)    = [];
infraPairCorrMidR(~goodRunsR|singleChRow,:)      = [];
infraPairCorrDeepR(~goodRunsR|singleChRow,:)     = [];

medPairASuperBMidR(~goodRunsR|singleChRow,:)  = [];
medPairASuperBDeepR(~goodRunsR|singleChRow,:) = [];
medPairAMidBSuperR(~goodRunsR|singleChRow,:)  = [];
medPairAMidBDeepR(~goodRunsR|singleChRow,:)   = [];
medPairADeepBSuperR(~goodRunsR|singleChRow,:) = [];
medPairADeepBMidR(~goodRunsR|singleChRow,:)   = [];

envelopeASuperBMidR(~goodRunsR|singleChRow,:)  = [];
envelopeASuperBDeepR(~goodRunsR|singleChRow,:) = [];
envelopeAMidBSuperR(~goodRunsR|singleChRow,:)  = [];
envelopeAMidBDeepR(~goodRunsR|singleChRow,:)   = [];
envelopeADeepBSuperR(~goodRunsR|singleChRow,:) = [];
envelopeADeepBMidR(~goodRunsR|singleChRow,:)   = [];

infraASuperBMidR(~goodRunsR|singleChRow,:)  = [];
infraASuperBDeepR(~goodRunsR|singleChRow,:) = [];
infraAMidBSuperR(~goodRunsR|singleChRow,:)  = [];
infraAMidBDeepR(~goodRunsR|singleChRow,:)   = [];
infraADeepBSuperR(~goodRunsR|singleChRow,:) = [];
infraADeepBMidR(~goodRunsR|singleChRow,:)   = [];

intraCorrAR(~goodRunsR|singleChRow,:)        = [];
intraCorrBR(~goodRunsR|singleChRow,:)        = [];
envelopeIntraAAllR(~goodRunsR|singleChRow,:) = [];
envelopeIntraBAllR(~goodRunsR|singleChRow,:) = [];
infraIntraAAllR(~goodRunsR|singleChRow,:)    = [];
infraIntraBAllR(~goodRunsR|singleChRow,:)    = [];

%% Get cross-frequency correlations
fs        = 1e3; % Sampling frequency
thetaBand = [6 8];    [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters
alphaBand = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass');
betaBand  = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass');

% Check if certain variables are stored....
if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') 
    varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    varIdx  = ismember('medPairCorrFreq',varInfo);
else 
    varIdx = 1;
end
% Get the combination
freqCombs = nchoosek(1:size(bandLabels,2),2);

% Calculate pairwise correlations between and within probes
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') || ~varIdx
    tic; 
    for iDate = 1: size(allDates,1)
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

            rawA   = allProbeData{fileNum,iDate}.raw1Ch;
            rawB   = allProbeData{fileNum,iDate}.raw2Ch;

            % Remove bad channels
            probeA(:,badElecA{fileNum,iDate}) = [];
            probeB(:,badElecB{fileNum,iDate}) = [];

            rawA(:,badElecA{fileNum,iDate}) = [];
            rawB(:,badElecB{fileNum,iDate}) = [];

            % Remove bad times
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];

            rawA(allBadTimes{fileNum,iDate},:) = [];
            rawB(allBadTimes{fileNum,iDate},:) = [];

            chA = estChInCortexA{iDate}(iRun,:);
            chB = estChInCortexB{iDate}(iRun,:);

            if chA(1) == 0 || chB(1) == 0 || isempty(probeA) || isempty(probeB)
                medPairCorrFreqTime(1:size(freqCombs,2),iDate,iRun)= NaN;
                medPairCorrFreqPower(1:size(freqCombs,2),iDate,iRun)= NaN;
                medPairCorrFreqInfra(1:size(freqCombs,2),iDate,iRun)= NaN;

                % medPACFreqTime(1:size(freqCombs,2),iDate,iRun)= NaN;
                % medPACFreqPower(1:size(freqCombs,2),iDate,iRun)= NaN;
                % medPACFreqInfra(1:size(freqCombs,2),iDate,iRun)= NaN;
                
                continue;
            end
          
            for iComb = 1:size(freqCombs,1)
                clear combVals pAFreq pBFreq envelopeABand envelopeBBand
                combVals = freqCombs(iComb,:);

                for iPair = 1:length(combVals)
                    switch combVals(iPair)
                        case 1
                            pAFreq(:,:,iPair) = single(filtfilt(bT,aT,double(probeA(:,chA(1):chA(2)))));
                            pBFreq(:,:,iPair) = single(filtfilt(bT,aT,double(probeB(:,chB(1):chB(2)))));
                        case 2
                            pAFreq(:,:,iPair) = single(filtfilt(bA,aA,double(probeA(:,chA(1):chA(2)))));
                            pBFreq(:,:,iPair) = single(filtfilt(bA,aA,double(probeB(:,chB(1):chB(2)))));
                        case 3
                            pAFreq(:,:,iPair) = single(filtfilt(bB,aB,double(probeA(:,chA(1):chA(2)))));
                            pBFreq(:,:,iPair) = single(filtfilt(bB,aB,double(probeB(:,chB(1):chB(2)))));
                        case 4
                            pAFreq(:,:,iPair) = single(filtfilt(bG,aG,double(probeA(:,chA(1):chA(2)))));
                            pBFreq(:,:,iPair) = single(filtfilt(bG,aG,double(probeB(:,chB(1):chB(2)))));
                        case 5
                            pAFreq(:,:,iPair) = rawA(:,chA(1):chA(2));
                            pBFreq(:,:,iPair) = rawB(:,chB(1):chB(2));
                    end

                    % Get instantaneous power and correlate the powers
                    envelopeABand(:,:,iPair) = envelope(abs(pAFreq(:,:,iPair)),5);
                    envelopeBBand(:,:,iPair) = envelope(abs(pBFreq(:,:,iPair)),5);
                end

                medPairCorrFreqTime(iComb,iDate,iRun) = median([corr(squeeze(pAFreq(:,:,1)),squeeze(pBFreq(:,:,2))) corr(squeeze(pAFreq(:,:,2)),squeeze(pBFreq(:,:,1)))],'all','omitnan');
              

                % Correlate instantaneous band power
                medPairCorrFreqPower(iComb,iDate,iRun) = median([corr(squeeze(envelopeABand(:,:,1)),squeeze(envelopeBBand(:,:,2))) corr(squeeze(envelopeABand(:,:,2)),squeeze(envelopeBBand(:,:,1)))],'all','omitnan');

                enSizeA = size(envelopeABand);
                enSizeB = size(envelopeBBand);

                % Filter envelope from 0.01-0.1 Hz
                envelopeABand = [envelopeABand; envelopeABand ;envelopeABand];
                infraSlowA = filtfilt(sos,g,double(envelopeABand));
                infraSlowA = single(infraSlowA(enSizeA(1)+1:(end-enSizeA(1)),:,:));

                envelopeBBand = [envelopeBBand; envelopeBBand ;envelopeBBand];
                infraSlowB = filtfilt(sos,g,double(envelopeBBand));
                infraSlowB = single(infraSlowB(enSizeB(1)+1:(end-enSizeB(1)),:,:));

                % Correlate infraslow flucutuations in instantaneous band power
                medPairCorrFreqInfra(iComb,iDate,iRun) =  median([corr(squeeze(infraSlowA(:,:,1)),squeeze(infraSlowB(:,:,2))) corr(squeeze(infraSlowA(:,:,2)),squeeze(infraSlowB(:,:,1)))],'all','omitnan');
            end
        end
    end
    
    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'medPairCorrFreqTime','-append');

toc;
else
    disp('Loading saved variables...')
    clear allVars
    allVars = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    fieldNames = fieldnames(allVars); 

    for iL = 1:length(fieldNames)
        nameVal = fieldNames{iL};
        assignin('base',nameVal,allVars.(nameVal));
    end
    clear allVars;
end

medPairCorrFreqTimeR = reshape(medPairCorrFreqTime,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]); 
medPairCorrFreqTimeR(:,removeDataIdx)          = [];
medPairCorrFreqTimeR(:,~goodRunsR) = [];

medPairCorrFreqPowerR = reshape(medPairCorrFreqPower,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]); 
medPairCorrFreqPowerR(:,removeDataIdx)          = [];
medPairCorrFreqPowerR(:,~goodRunsR) = [];

medPairCorrFreqInfraR = reshape(medPairCorrFreqInfra,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]); 
medPairCorrFreqInfraR(:,removeDataIdx)          = [];
medPairCorrFreqInfraR(:,~goodRunsR) = [];

freqCombNames = {'Theta-Alpha','Theta-Beta','Theta-Gamma','Theta-Spiking',...
    'Alpha-Beta','Alpha-Gamma','Alpha-Spiking','Beta-Gamma','Beta-Spiking','Gamma-Spiking'}; 
%%
for iPlot = 1:3
    switch iPlot 
        case 1
            plotVar = medPairCorrFreqTimeR;
            plotTitle = 'Time series';
        case 2
            plotVar = medPairCorrFreqPowerR;
            plotTitle = 'Power';
        case 3
            plotVar = medPairCorrFreqInfraR; 
            plotTitle = 'Infraslow'; 
    end
    % figure; xVal = connValsR; 
    % for iFig = 1: size(medPairCorrFreqTimeR,1)
    %     subplot(2,5,iFig);
    %     yVal = plotVar(iFig,:);
    %     plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
    %     coeff = polyfit(xVal,yVal,1);
    %     xFit  = linspace(min(xVal),max(xVal),1000);
    %     yFit  = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
    %     plot(xFit,yFit,'-k','LineWidth',1);
    %     text(0.7, 0.35,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    %     text(0.7, 0.3,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
    %     title(freqCombNames{iFig}); ylim([-0.4 0.4]); xlim([-0.6 1]); axis square;
    %     if iFig== 1; xlabel('Functional connectivity'); ylabel('Cross-frequency correlations'); end 
    % end
    % sgtitle(plotTitle);

     figure; xVal = distValsR; 
    for iFig = 1: size(medPairCorrFreqTimeR,1)
        subplot(2,5,iFig);
        yVal = plotVar(iFig,:)';
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
        text(10, 0.35,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
        text(10, 0.3,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
        title(freqCombNames{iFig}); ylim([-0.4 0.4]); xlim([0 20]); axis square;
        if iFig== 1; xlabel('Distance (mm)'); ylabel('Cross-frequency correlations'); end
    end
    sgtitle(plotTitle);
end



%% Cross-frequency coupling - Phase Amplitude coupling between frequencies
for iDate = 1: size(allDates,1)
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
        if chA(1)== 0 || chB(1)==0; continue; end

        pAGamma = single(filtfilt(bG,aG,double(probeA(:,chA(1):chA(2)))));
        pBGamma = single(filtfilt(bG,aG,double(probeB(:,chB(1):chB(2)))));

        pATheta = single(filtfilt(bT,aT,double(probeA(:,chA(1):chA(2)))));
        pBTheta = single(filtfilt(bT,aT,double(probeB(:,chB(1):chB(2)))));
        
        pAAlpha = single(filtfilt(bA,aA,double(probeA(:,chA(1):chA(2)))));
        pBAlpha = single(filtfilt(bA,aA,double(probeB(:,chB(1):chB(2)))));

        pABeta = single(filtfilt(bB,aB,double(probeA(:,chA(1):chA(2)))));
        pBBeta = single(filtfilt(bB,aB,double(probeB(:,chB(1):chB(2)))));
      
        for iFreq = 1:3
            clear combVals pAFreq pBFreq envelopeABand envelopeBBand

            switch iFreq
                case 1
                    pAVal = pATheta;
                    pBVal = pBTheta;
                case 2
                    pAVal = pAAlpha;
                    pBVal = pBAlpha;
                case 3
                    pAVal = pABeta;
                    pBVal = pBBeta;
            end

            numPoints = size(pAGamma,1);
            numSurrogate = 200; % to get the joint distribution
            minSkip  = fs;
            maxSkip  = numPoints-fs;
            skip     = ceil(numPoints.*rand(numSurrogate*2,1));
            skip(skip>maxSkip) = [];
            skip(skip<minSkip) = [];

            for iRef = 1:2
                switch iRef
                    case 1
                        pRef = pAGamma; 
                        pMov = pBVal;
                    case 2
                        pRef = pBGamma; 
                        pMov = pAVal; 
                end
                clear amplitude phase z surrAmplitude surrogate surrMean surrStd
                numElec   = min([size(pRef,2) size(pMov,2)]);
                surrogate = zeros(numSurrogate,numElec);

                % Get analytical signal
                amplitude = abs(hilbert(pAGamma(:,1:numElec))); % Gamma
                phase     = angle(hilbert(pBTheta(:,1:numElec))); % Theta

                % Complex valued signal
                z = amplitude.*exp(1i*phase);

                % Get joint distribution and calculate mean of z to get modulation
                % index
                mRaw = mean(z);
                for iS = 1:numSurrogate
                    surrAmplitude = [amplitude(skip(iS):end,:); amplitude(1:skip(iS)-1,:)];
                    surrogate(iS,:) = abs(mean(surrAmplitude.*exp(1i*phase)));
                end

                % Fit Gaussian to surrogate data
                [surrMean, surrStd] = normfit(surrogate);
                mNormLen(iDate,iRun,iRef,iFreq) = median(abs((mRaw)-surrMean)./surrStd);
                % mNormPhase = angle(mRaw);
                % mNorm = mNormLen*exp(1i*mNormPhase);
            end
        end
    end
end
mNormNew = reshape(squeeze(mNormLen(:,:,1,:)),[40 3]);
connValsNew = connValsR(1:40);
isZero = mNormNew(:,1)== 0;
mNormNew(isZero,:) = [];
connValsNew(isZero) = [];

figure; scatter(connValsNew,mNormNew); 


%% Plotting
% Pairwise correlations vs functional connectivity
showScatterPlots(connValsR,medPairCorrR,medCorrEnvelopeR,medCorrInfraSlowR,...
    'Functional Connectivity','Pairwise correlations',[-0.6 1],[-0.6 1],...
    0.8,-0.5,-0.4,bandLabels);

%% Pairwise correlations vs distance
showScatterPlots(distValsR,medPairCorrR,medCorrEnvelopeR,medCorrInfraSlowR,'Distance',...
    'Pairwise correlations',[0 20],[-0.6 1],12,-0.5,-0.4,bandLabels);

%% Distributions of pairwise correlations and distance
figure; subplot(121); histogram(connValsR,-0.6:0.2:1); box off;xlabel('Functional connectivity');ylim([0 30]);xticks(-0.6:0.2:1);
subplot(122); histogram(distValsR,0:2:16); xlabel('Distance (mm)'); box off;ylim([0 30]); xticks(0:2:20);

%% Partial correlation calculations
for iType = 1:3
    clear var
    switch iType
        case 1
            var = medPairCorrR;
        case 2
            var = medCorrEnvelopeR;
        case 3
            var = medCorrInfraSlowR; 
    end
    for iBand = 1:5
        [rhoVal(iType,iBand),pVal(iType,iBand)] = partialcorr(connValsR,var(:,iBand),distValsR);
    end
end

%%  Show the distributions for within probe pairwise correlations
oneIdx = find(medIntraCorrBR(:,1)==1);
medIntraCorrBR(oneIdx,:) = [];
envelopeIntraCorrBR(oneIdx,:) = [];
infraIntraCorrBR(oneIdx,:) = [];

figure;
plotIdxConn = 1; 
for iType = 1:3
    switch iType
        case 1
            plotVal   = [medIntraCorrAR; medIntraCorrBR];
            plotLabel = 'Timeseries';
        case 2
            plotVal   = [envelopeIntraCorrAR; envelopeIntraCorrBR];
            plotLabel = 'Envelope';
        case 3
            plotVal   = [infraIntraCorrAR; infraIntraCorrBR];
            plotLabel = 'Infraslow';
    end
    for iBand = 1:5
        subplot(3,5,plotIdxConn);
        h = histogram(plotVal(:,iBand),0:0.1:1); 
        xlabel('Within probe correlations'); ylabel('Count');
        xticks(0:0.2:1);

        if plotIdxConn == 5 || plotIdxConn == 10
            ylim([0 100]); yticks(0:20:150);
        % elseif 
        %     ylim([0 100]); yticks(0:10:100);
        else
             ylim([0 40]); yticks(0:10:60);
        end      
        title([plotLabel ' - ' bandLabels{iBand}]); box off;
        plotIdxConn = plotIdxConn+1;
    end
end
%% Condensing all relationships into scatters

for iType = 1:3
    switch iType
        case 1
            vals = medPairCorrR;
        case 2
            vals = medCorrEnvelopeR;
        case 3
            vals = medCorrInfraSlowR;
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
    y = zscore(connValsR);
    X  = zscore(vals);
    clear  beta sigma eVal covMat

    % Get the linear model 
    mdl{iType}        = fitlm(X,y,'VarNames',{'Theta','Alpha','Beta','Gamma','Spk','FC'});
    
    % Relative weights analysis to determine contributions of predictors
    [relImp(iType), r2(iType)] = rwa(X,y,bandLabels');

    % Dominance analysis
    [relativeImportance(:,iType),rsqDominance(iType)] = dominance(X,y);
    percentImportance(:,iType) = 100*relativeImportance(:,iType)./rsqDominance(iType);
end

% Plot the dominance weights
figure;bar(percentImportance);xticks(1:5);xticklabels(bandLabels); ylabel('Relative importance to overall R^2 (%)');
legend(timeLabels,'Location','northeast'); box off; title('Dominance weights');

% Plot the relative weights
figure; bar([relImp.scores]);xticks(1:5);xticklabels(bandLabels); ylabel('Relative weights to overall R^2 (%)');
legend(timeLabels,'Location','northeast'); box off; title('Relative weights'); 

% Plot the beta values obtained from the Relative weights method
figure; bar([relImp.betaVals]);xticklabels(bandLabels); ylabel('Standardized beta values'); 
legend(timeLabels,'Location','northeast'); box off; title('Relative weights');colororder('meadow')

% Bootstrapping relative weights and/or dominance analysis
tic;
for iType = 1:3
    clear allVar X y 
    switch iType
        case 1
            vals = medPairCorrR;
        case 2
            vals = medCorrEnvelopeR;
        case 3
            vals = medCorrInfraSlowR;
    end
    y = zscore(connValsR);
    X  = zscore(vals);
    allVar  = [X y];

    % Perform resampling and get relative weights or dominance weights
    rowIdx  = 1:size(allVar,1); 
    numRows = length(rowIdx);
  
    for iRep = 1:100
        clear newRows newData  
        newRows = datasample(rowIdx,numRows,'Replace',true);
        newData = allVar(newRows,:);
        
        % Relative weights analysis 
        [relImpBoot(iRep,iType), r2Boot(iRep,iType)] = rwa(newData(:,1:5),newData(:,6),bandLabels');

        % Dominance analysis
        [domWeightsBoot(:,iRep,iType),rsqDominanceBoot(iRep,iType)] = dominance(newData(:,1:5),newData(:,6));
        percentImpBoot(:,iRep,iType) = 100*domWeightsBoot(:,iRep,iType)./rsqDominanceBoot(iRep,iType);
    end   
end
toc;

% Dominance statistics 
meanPercentImpBoot = squeeze(mean(percentImpBoot,2,'omitnan')); 
stdErrorPerImpBoot = squeeze(std(percentImpBoot,[],2,'omitnan'))./sqrt(size(percentImpBoot,2));
r2BootMean         = mean(r2Boot,1,'omitnan');

% Plot the dominance weights
figure;b = bar(meanPercentImpBoot);xticks(1:5);xticklabels(bandLabels); ylabel('Relative importance to overall R^2 (%)');
legend(timeLabels,'Location','northeast','AutoUpdate','off'); box off; title('Bootstrapped dominance weights'); hold on; 
for iB = 1:numel(b)
    xtips = b(iB).XEndPoints;
    ytips = b(iB).YEndPoints; 
    errorbar(xtips,ytips,2.*stdErrorPerImpBoot(:,iB),'.k','Markersize',0.1);
end

% RWA stats
meanRelImp     = [mean([relImpBoot(:,1).scores],2) mean([relImpBoot(:,2).scores],2) mean([relImpBoot(:,3).scores],2) ]; 
stdErrorRelImp = [std([relImpBoot(:,1).scores],[],2) std([relImpBoot(:,2).scores],[],2) std([relImpBoot(:,3).scores],[],2)]./sqrt(size(percentImpBoot,2));

figure; b= bar(meanRelImp);xticks(1:5);xticklabels(bandLabels); ylabel('Relative weights to overall R^2 (%)');
legend(timeLabels,'Location','northeast','AutoUpdate','off'); box off; title('Bootstrapped Relative weights'); hold on; 
for iB = 1:numel(b)
    xtips = b(iB).XEndPoints;
    ytips = b(iB).YEndPoints; 
    errorbar(xtips,ytips,2.*stdErrorRelImp(:,iB),'.k','Markersize',0.1);
end

meanRelWgts    = [mean([relImpBoot(:,1).betaVals],2) mean([relImpBoot(:,2).betaVals],2) mean([relImpBoot(:,3).betaVals],2) ]; 
stdErrorRelWgt = [std([relImpBoot(:,1).betaVals],[],2) std([relImpBoot(:,2).betaVals],[],2) std([relImpBoot(:,3).betaVals],[],2)]./sqrt(size(percentImpBoot,2));

% Plot the beta values obtained from the Relative weights method
figure;b= bar(meanRelWgts);xticklabels(bandLabels); ylabel('Standardized beta values'); hold on; 
legend(timeLabels,'Location','northeast','AutoUpdate','off'); box off; title('Bootstrapped Relative weights');colororder('meadow')
for iB = 1:numel(b)
    xtips = b(iB).XEndPoints;
    ytips = b(iB).YEndPoints; 
    errorbar(xtips,ytips,2.*stdErrorRelWgt(:,iB),'.k','Markersize',0.1);
end


% %% Partial regression
%  figure;
%  clear xl yl xs ys beta pctvar mse statsPar vipScore
% for iType = 1:3
%     switch iType
%         case 1
%                vals = medPairCorrR;
%         case 2
%                vals = medCorrEnvelopeR;
%         case 3
%                vals = medCorrInfraSlowR;
%     end
%     for iPred = 1:size(vals,2)
%         otherIdx            = setdiff(1:size(vals,2),iPred);
%         mdlYOthers{iPred,iType}   = fitlm(vals(:,otherIdx),connValsR);
%         residualsY(:,iPred,iType) = mdlYOthers{iPred,iType}.Residuals.Raw;
% 
%         mdlXOthers{iPred,iType}   = fitlm(vals(:,otherIdx),vals(:,iPred));
%         residualsX(:,iPred,iType) = mdlXOthers{iPred,iType}.Residuals.Raw;
%         tempMdl    = fitlm(residualsX(:,iPred,iType),residualsY(:,iPred,iType));
%         slopeMdl(iPred,iType) = tempMdl.Coefficients.Estimate(2);
%     end
% 
%     [xl,yl,xs,ys,beta{iType},pctvar(:,:,iType),mse{iType},...
%         statsPar{iType}] = plsregress(vals,connValsR,5);
% 
%     w0 = statsPar{iType}.W ./ sqrt(sum(statsPar{iType}.W.^2,1));
% 
%     p = size(xl,1);
%     sumSq = sum(xs.^2,1).*sum(yl.^2,1);
%     vipScore(:,iType) = sqrt(p* sum(sumSq.*(w0.^2),2) ./ sum(sumSq,2));
% 
%         % indVIP = find(vipScore >= 1);
% 
%     % Figures
% 
%     subplot(3,5,1+5*(iType-1)); scatter(residualsY(:,1,iType),residualsX(:,1,iType),'filled'); ylim([-0.4 0.4]); xlim([-0.4 0.4]);
%     subplot(3,5,2+5*(iType-1)); scatter(residualsY(:,2,iType),residualsX(:,2,iType),'filled'); ylim([-0.4 0.4]); xlim([-0.4 0.4]);
%     subplot(3,5,3+5*(iType-1)); scatter(residualsY(:,3,iType),residualsX(:,3,iType),'filled'); ylim([-0.4 0.4]); xlim([-0.4 0.4]);
%     subplot(3,5,4+5*(iType-1)); scatter(residualsY(:,4,iType),residualsX(:,4,iType),'filled'); ylim([-0.4 0.4]); xlim([-0.4 0.4]);
%     subplot(3,5,5+5*(iType-1)); scatter(residualsY(:,5,iType),residualsX(:,5,iType),'filled'); ylim([-0.4 0.4]); xlim([-0.4 0.4]);
% end
% 
% figure; 
% plot(1:5,cumsum(100*squeeze(pctvar(2,:,1))),'-'); hold on;
% plot(1:5,cumsum(100*squeeze(pctvar(2,:,2))),'-'); hold on;
% plot(1:5,cumsum(100*squeeze(pctvar(2,:,3))),'-'); box off;
% legend('Time-series','Power','Infraslow power','location','northwest');
% xlabel('Components/Predictors'); ylabel('% variance explained'); xticks(1:5); 
% 
% % VIP score
% figure;
% plot(1:5,vipScore,'.','MarkerSize',20);
% legend('Time-series','Power','Infraslow power','location','northwest','autoupdate','off');
% xlabel('Components/Predictors'); ylabel('VIP score'); xticks(1:5); box off;
% yline(1,'Linewidth',1);


%% Cross-laminar correlations vs Functional connectivity
showScatterPlotsLayers(connValsR,medPairCorrSuperR,medPairCorrMidR,medPairCorrDeepR,...
    envelopePairCorrSuperR,envelopePairCorrMidR,envelopePairCorrDeepR,...
    infraPairCorrSuperR,infraPairCorrMidR,infraPairCorrDeepR,...
    'Functional Connectivity','Pairwise correlations',[-0.6 1],[-0.6 1],...
    0.8,-0.5,-0.4,bandLabels);

%% Cross-laminar correlations vs Distance
showScatterPlotsLayers(distValsR,medPairCorrSuperR,medPairCorrMidR,medPairCorrDeepR,...
    envelopePairCorrSuperR,envelopePairCorrMidR,envelopePairCorrDeepR,...
    infraPairCorrSuperR,infraPairCorrMidR,infraPairCorrDeepR,...
    'Distance','Pairwise correlations',[0 20],[-0.6 1],12,-0.5,-0.4,bandLabels);


%% Plot the between compartment correlations
connValsNew = connValsR;
connValsNew(singleChRow) = [];

for iPlot = 1: 3
    clear yValsSuperMid yValsSuperDeep yValsMidDeep 
    switch iPlot
        case 1
            yValsSuperMid   = (medPairASuperBMidR + medPairAMidBSuperR)./2;
            yValsSuperDeep  = (medPairASuperBDeepR+ medPairADeepBSuperR)./2;
            yValsMidDeep    = (medPairAMidBDeepR+ medPairADeepBMidR)./2;
            yValsSuperSuper = medPairCorrSuperR;
            yValsDeepDeep   = medPairCorrDeepR;
            yValsMidMid     = medPairCorrMidR;
            typeVal         = 'Time series';

        case 2
            yValsSuperMid   = (envelopeASuperBMidR+ envelopeAMidBSuperR)./2;
            yValsSuperDeep  = (envelopeASuperBDeepR+ envelopeADeepBSuperR)./2;
            yValsMidDeep    = (envelopeAMidBDeepR+envelopeADeepBMidR)./2;
            yValsSuperSuper = envelopePairCorrSuperR;
            yValsDeepDeep   = envelopePairCorrDeepR;
            yValsMidMid     = envelopePairCorrMidR;
            typeVal         = 'Envelope';

        case 3
            yValsSuperMid   = (infraASuperBMidR+ infraAMidBSuperR)./2;
            yValsSuperDeep  = (infraASuperBDeepR+ infraADeepBSuperR)./2;
            yValsMidDeep    = (infraAMidBDeepR+ infraADeepBMidR)./2;
            yValsSuperSuper = infraPairCorrSuperR;
            yValsDeepDeep   = infraPairCorrDeepR;
            yValsMidMid     = infraPairCorrMidR;
            typeVal         = 'Infraslow';
    end

    % figure;%('position',[932,103,709,838]);
    corrSuperMid(:,iPlot)  = corr(yValsSuperMid,connValsNew);
    corrSuperDeep(:,iPlot) = corr(yValsSuperDeep,connValsNew);
    corrMidDeep(:,iPlot) = corr(yValsMidDeep,connValsNew);

    corrSuperSuper(:,iPlot) = corr(yValsSuperSuper,connValsNew);
    corrDeepDeep(:,iPlot) = corr(yValsDeepDeep,connValsNew);
    corrMidMid(:,iPlot) = corr(yValsMidMid,connValsNew);

    % for iBand = 1:5 % Feedback
    %     subplot(3,5,idx); showLinearFit(connValsNew,yValsSuperDeep(:,iBand),0.8,-0.5,-0.4);
    %     title([bandLabels{iBand} ' - Super x Deep']);xlim([-0.6 1]); ylim([-0.6 1]);axis square;
    % 
    %     subplot(3,5,idx+5); showLinearFit(connValsR,yValsSuperSuper(:,iBand),0.8,-0.5,-0.4);
    %     title([bandLabels{iBand} ' - Super x Super ']); xlim([-0.6 1]); ylim([-0.6 1]); axis square;
    % 
    %     subplot(3,5,idx+10); showLinearFit(connValsR,yValsDeepDeep(:,iBand),0.8,-0.5,-0.4);
    %     title([bandLabels{iBand} ' - Deep  x Deep ']);xlim([-0.6 1]); ylim([-0.6 1]);axis square;
    % 
    %     if iBand==1; xlabel('Functional connectivity'); ylabel('Pairwise correlations');end
    %     sgtitle(['Feedback: ' typeVal]);
    % 
    %     % [rhoValSM(iPlot,iBand),pValSM(iType,iBand)] = partialcorr(connValsNew,yValsSuperMid(:,iBand),distValsR);
    %     % [rhoValSD(iPlot,iBand),pValSD(iType,iBand)] = partialcorr(connValsNew,yValsSuperDeep(:,iBand),distValsR);
    %     % [rhoValMD(iPlot,iBand),pValMD(iType,iBand)] = partialcorr(connValsNew,yValsMidDeep(:,iBand),distValsR);
    %     idx = idx+1;
    %     % if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs'],'dir')
    %     %     [~,~]= mkdir(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs']);
    %     % end
    %     % f = gcf;
    %     % exportgraphics(f,['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs\' bandLabels{iBand} '-' typeVal '.png'],'Resolution',300);
    %     % close gcf;
    % 
    % end
end
timeLabels = {'Time series','Power','Infraslow'};
figure;
subplot(321); imagesc(imgaussfilt(corrSuperMid));   title('Super x Middle');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(322); imagesc(imgaussfilt(corrSuperDeep));  title('Super x Deep');    yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(323); imagesc(imgaussfilt(corrMidDeep));    title('Deep x Middle');   yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(324); imagesc(imgaussfilt(corrSuperSuper)); title('Super x Super');   yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(325); imagesc(imgaussfilt(corrDeepDeep));   title('Middle x Middle'); yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(326); imagesc(imgaussfilt(corrMidMid));     title('Deep x Deep');     yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;

%% Plot the distributions of pairwise correlations  
plotIdx =1; figure('units','normalized','outerposition',[0 0 1 1]);
for iPlot = 1: 3
    switch iPlot
        case 1
            yValsSuperMid   = [medPairASuperBMidR ;medPairAMidBSuperR];
            yValsSuperDeep  = [medPairASuperBDeepR; medPairADeepBSuperR];
            yValsMidDeep    = [medPairAMidBDeepR; medPairADeepBMidR];
            yValsSuperSuper = [medPairCorrSuperR; NaN(size(medPairCorrSuperR))];
            yValsDeepDeep   = [medPairCorrDeepR;NaN(size(medPairCorrSuperR))];
            yValsMidMid     = [medPairCorrMidR;NaN(size(medPairCorrSuperR))];
            typeVal         = 'Time series';

        case 2
            yValsSuperMid   = [envelopeASuperBMidR; envelopeAMidBSuperR];
            yValsSuperDeep  = [envelopeASuperBDeepR; envelopeADeepBSuperR];
            yValsMidDeep    = [envelopeAMidBDeepR;  envelopeADeepBMidR];
            yValsSuperSuper = [envelopePairCorrSuperR;NaN(size(medPairCorrSuperR))];
            yValsDeepDeep   = [envelopePairCorrDeepR;NaN(size(medPairCorrSuperR))];
            yValsMidMid     = [envelopePairCorrMidR;NaN(size(medPairCorrSuperR))];
            typeVal         = 'Envelope';
        case 3
            yValsSuperMid   = [infraASuperBMidR; infraAMidBSuperR];
            yValsSuperDeep  = [infraASuperBDeepR; infraADeepBSuperR];
            yValsMidDeep    = [infraAMidBDeepR;  infraADeepBMidR];
            yValsSuperSuper = [infraPairCorrSuperR;NaN(size(medPairCorrSuperR))];
            yValsDeepDeep   = [infraPairCorrDeepR;NaN(size(medPairCorrSuperR))];
            yValsMidMid     = [infraPairCorrMidR;NaN(size(medPairCorrSuperR))];
            typeVal         = 'Infraslow';
    end

    for iBand = 1:5
        subplot(3,5,plotIdx); 
        boxplot([yValsSuperMid(:,iBand) yValsSuperDeep(:,iBand) yValsMidDeep(:,iBand)...
            yValsSuperSuper(:,iBand) yValsMidMid(:,iBand) yValsDeepDeep(:,iBand)],{'Super x Middle',...
            'Super x Deep','Middle x Deep','Super x Super','Mid x Mid', 'Deep x Deep'});
        title([typeVal '- ' bandLabels{iBand}]); box off; ylim([-0.6 1]); yticks(-0.6:0.2:1);%axis square;
        plotIdx = plotIdx+1;
    end
end

%% Separating the correlations based on the "reference" probe

for iPlot = 1:3
    switch iPlot
        case 1
            corrASuperBMid(:,iPlot)  = corr(medPairASuperBMidR,connValsNew);
            corrBSuperAMid(:,iPlot)  = corr(medPairAMidBSuperR,connValsNew);
            corrASuperBDeep(:,iPlot) = corr(medPairASuperBDeepR,connValsNew);
            corrBSuperADeep(:,iPlot) = corr(medPairADeepBSuperR,connValsNew);
            corrAMidBDeep(:,iPlot)   = corr(medPairAMidBDeepR,connValsNew);
            corrBMidADeep(:,iPlot)   = corr(medPairADeepBMidR,connValsNew);
            typeVal                  = 'Time series';

        case 2
            corrASuperBMid(:,iPlot)  = corr(envelopeASuperBMidR,connValsNew);
            corrBSuperAMid(:,iPlot)  = corr(envelopeAMidBSuperR,connValsNew);
            corrASuperBDeep(:,iPlot) = corr(envelopeASuperBDeepR,connValsNew);
            corrBSuperADeep(:,iPlot) = corr(envelopeADeepBSuperR,connValsNew);
            corrAMidBDeep(:,iPlot)   = corr(envelopeAMidBDeepR,connValsNew);
            corrBMidADeep            = corr(envelopeADeepBMidR,connValsNew);
            typeVal                  = 'Envelope';
 
        case 3
            corrASuperBMid(:,iPlot)  = corr(infraASuperBMidR,connValsNew);
            corrBSuperAMid(:,iPlot)  = corr(infraAMidBSuperR,connValsNew);
            corrASuperBDeep(:,iPlot) = corr(infraASuperBDeepR,connValsNew);
            corrBSuperADeep(:,iPlot) = corr(infraADeepBSuperR,connValsNew);
            corrAMidBDeep(:,iPlot)   = corr(infraAMidBDeepR,connValsNew);
            corrBMidADeep(:,iPlot)   = corr(infraADeepBMidR,connValsNew);
            typeVal                  = 'Infraslow';
    end
end

figure;
subplot(3,2,1); imagesc(imgaussfilt(corrASuperBMid));   title('A Super x B Middle');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(3,2,2); imagesc(imgaussfilt(corrBSuperAMid));   title('B Super x A Middle');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(3,2,3); imagesc(imgaussfilt(corrASuperBDeep));   title('A Super x B Deep');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(3,2,4); imagesc(imgaussfilt(corrBSuperADeep));   title('B Super x A Deep');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(3,2,5); imagesc(imgaussfilt(corrAMidBDeep));   title('A Middle x B Deep');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;
subplot(3,2,6); imagesc(imgaussfilt(corrBMidADeep));   title('B Middle x A Deep');  yticklabels(bandLabels); xticks(1:3); xticklabels(timeLabels); clim([-0.2 0.5]); colormap jet; axis square; colorbar;

%% Plot correlations between compartments in the same probe
chSplit = 6; 
superCh = 1:chSplit; 
midCh   = chSplit+1:chSplit*2;

% Remove recordings with only a single channel
% singleChRow = cellfun(@(x) isscalar(x),intraCorrBR(:,1));
intraCorrAR(singleChRow,:) = [];
intraCorrBR(singleChRow,:) = [];

envelopeIntraAAllR(singleChRow,:) = [];
envelopeIntraBAllR(singleChRow,:) = [];

infraIntraAAllR(singleChRow,:) = [];
infraIntraBAllR(singleChRow,:) = [];

% Get median correlations between compartments
superMidTimeSeriesA = cell2mat(cellfun(@(x) median(x(midCh,superCh),'all','omitnan'),intraCorrAR,'un',0));
superMidTimeSeriesB = cell2mat(cellfun(@(x) median(x(midCh,superCh),'all','omitnan'),intraCorrBR,'un',0));

superDeepTimeSeriesA = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,superCh),'all','omitnan'),intraCorrAR,'un',0));
superDeepTimeSeriesB = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,superCh),'all','omitnan'),intraCorrBR,'un',0));

midDeepTimeSeriesA = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,midCh),'all','omitnan'),intraCorrAR,'un',0));
midDeepTimeSeriesB = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,midCh),'all','omitnan'),intraCorrBR,'un',0));

superMidEnvelopeA = cell2mat(cellfun(@(x) median(x(midCh,superCh),'all','omitnan'),envelopeIntraAAllR,'un',0));
superMidEnvelopeB = cell2mat(cellfun(@(x) median(x(midCh,superCh),'all','omitnan'),envelopeIntraBAllR,'un',0));

superDeepEnvelopeA = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,superCh),'all','omitnan'),envelopeIntraAAllR,'un',0));
superDeepEnvelopeB = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,superCh),'all','omitnan'),envelopeIntraBAllR,'un',0));

midDeepEnvelopeA = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,midCh),'all','omitnan'),envelopeIntraAAllR,'un',0));
midDeepEnvelopeB = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,midCh),'all','omitnan'),envelopeIntraBAllR,'un',0));

superMidInfraA = cell2mat(cellfun(@(x) median(x(midCh,superCh),'all','omitnan'),infraIntraAAllR,'un',0));
superMidInfraB = cell2mat(cellfun(@(x) median(x(midCh,superCh),'all','omitnan'),infraIntraBAllR,'un',0));

superDeepInfraA = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,superCh),'all','omitnan'),infraIntraAAllR,'un',0));
superDeepInfraB = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,superCh),'all','omitnan'),infraIntraBAllR,'un',0));

midDeepInfraA = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,midCh),'all','omitnan'),infraIntraAAllR,'un',0));
midDeepInfraB = cell2mat(cellfun(@(x) median(x(chSplit*2+1:end,midCh),'all','omitnan'),infraIntraBAllR,'un',0));

%% Plotting per between compartment comparison across timescales
figure; t = tiledlayout(1,11);
for iFig = 1:3
    switch iFig
        case 1
            timeSeriesA = superMidTimeSeriesA;
            timeSeriesB = superMidTimeSeriesB;
            powerA      = superMidEnvelopeA;
            powerB      = superMidEnvelopeB;
            infraA      = superMidInfraA;
            infraB      = superMidInfraB;
            titleLabel  = 'Super x Middle';

        case 2
            timeSeriesA = superDeepTimeSeriesA;
            timeSeriesB = superDeepTimeSeriesB;
            powerA      = superDeepEnvelopeA;
            powerB      = superDeepEnvelopeB;
            infraA      = superDeepInfraA;
            infraB      = superDeepInfraB;
            titleLabel  = 'Super x Deep';

        case 3
            timeSeriesA = midDeepTimeSeriesA;
            timeSeriesB = midDeepTimeSeriesB;
            powerA      = midDeepEnvelopeA;
            powerB      = midDeepEnvelopeB;
            infraA      = midDeepInfraA;
            infraB      = midDeepInfraB;
            titleLabel  = 'Middle x Deep';
    end
    nexttile;
    imagesc([median(timeSeriesA,1,'omitnan'); median(timeSeriesB,1,'omitnan')]'); colormap jet; clim([-0.2 0.5]); 
    yticks(1:5); xticks(1:2); yticklabels(bandLabels);xticklabels({'A','B'}); axis square;
    title(['Timeseries - ' titleLabel]);

    nexttile;
    imagesc([median(powerA,1,'omitnan'); median(powerB,1,'omitnan')]'); colormap jet; clim([-0.2 0.5]);
    yticks(1:5); xticks(1:2); yticklabels([]);xticklabels({'A','B'});axis square;
    title('Power ');

    nexttile;
    imagesc([median(infraA,1,'omitnan'); median(infraB,1,'omitnan')]'); colormap jet; clim([-0.2 0.5]);
    yticks(1:5); xticks(1:2); yticklabels([]);xticklabels({'A','B'});axis square;
    title('Infraslow ');

    if iFig~=3
        nexttile;
        imagesc(); cla;
    else
        colorbar;
    end
end

delete(nexttile(4));
delete(nexttile(8));

%% Plotting each timescale separately
for iFig = 1:3
    switch iFig
        case 1
            superMidA  = superMidTimeSeriesA;
            superMidB  = superMidTimeSeriesB;
            superDeepA = superDeepTimeSeriesA;
            superDeepB = superDeepTimeSeriesB;
            midDeepA   = midDeepTimeSeriesA;
            midDeepB   = midDeepTimeSeriesB;
            typeLabel  = 'Time series';
      

        case 2
            superMidA  = superMidEnvelopeA;
            superMidB  = superMidEnvelopeB;
            superDeepA = superDeepEnvelopeA;
            superDeepB = superDeepEnvelopeB;
            midDeepA   = midDeepEnvelopeA;
            midDeepB   = midDeepEnvelopeB;
            typeLabel  = 'Power';

        case 3
            superMidA  = superMidInfraA;
            superMidB  = superMidInfraB;
            superDeepA = superDeepInfraA;
            superDeepB = superDeepInfraB;
            midDeepA   = midDeepInfraA;
            midDeepB   = midDeepInfraB;
            typeLabel  = 'Infraslow';
    end
    figure; 
    subplot(131); imagesc([median(superMidA,1,'omitnan'); median(superMidB,1,'omitnan')]'); colormap jet; clim([-0.2 0.5]);
    yticks(1:5); xticks(1:2); title('Super x Middle'); yticklabels(bandLabels);xticklabels({'A','B'});colorbar; axis square;

    subplot(132); imagesc([median(superDeepA,1,'omitnan'); median(superDeepB,1,'omitnan')]'); colormap jet; clim([-0.2 0.5]);
    yticks(1:5); xticks(1:2); title('Super x Deep'); yticklabels(bandLabels);xticklabels({'A','B'});colorbar; axis square;
        
    subplot(133); imagesc([median(midDeepA,1,'omitnan'); median(midDeepB,1,'omitnan')]'); colormap jet; clim([-0.2 0.5]);
    yticks(1:5); xticks(1:2); title('Middle x Deep'); yticklabels(bandLabels);xticklabels({'A','B'});colorbar; axis square;

    sgtitle(typeLabel);

    
    for iProbe = 1:2
        switch iProbe
            case 1
                superMid  = superMidA;
                superDeep = superDeepA;
                midDeep   = midDeepA;
                nameVal   = 'Probe A';
            case 2
                superMid  = superMidB;
                superDeep = superDeepB;
                midDeep   = midDeepB;
                nameVal   = 'Probe B';
        end

        figure; idx = 1;
        for iBand = 1:5
            subplot(5,3,idx); histogram(superMid(:,iBand),0:0.2:1); title(['Super x Middle: ' bandLabels{iBand}  ]); box off;
            subplot(5,3,idx+1); histogram(superDeep(:,iBand),0:0.2:1);title(['Super x Deep: ' bandLabels{iBand}]); box off;
            subplot(5,3,idx+2); histogram(midDeep(:,iBand),0:0.2:1);title(['Middle x Deep: ' bandLabels{iBand}]); box off;
            idx = idx+3; sgtitle([typeLabel ': ' nameVal]);
        end 
    end
end

%% Plot the between compartment correlations
distValsNew = [distValsR;distValsR];
for iPlot = 1: 3
    switch iPlot
        case 1
            yValsSuperMid   = [medPairASuperBMidR ;medPairAMidBSuperR];
            yValsSuperDeep  = [medPairASuperBDeepR; medPairADeepBSuperR];
            yValsMidDeep    = [medPairAMidBDeepR; medPairADeepBMidR];
            yValsSuperSuper = medPairCorrSuperR;
            yValsDeepDeep   = medPairCorrDeepR;
            typeVal         = 'Time series';

        case 2
            yValsSuperMid   = [envelopeASuperBMidR; envelopeAMidBSuperR];
            yValsSuperDeep  = [envelopeASuperBDeepR; envelopeADeepBSuperR];
            yValsMidDeep    = [envelopeAMidBDeepR;  envelopeADeepBMidR];
            yValsSuperSuper = envelopePairCorrSuperR;
            yValsDeepDeep   = envelopePairCorrDeepR;
            typeVal         = 'Envelope';
        case 3
            yValsSuperMid   = [infraASuperBMidR; infraAMidBSuperR];
            yValsSuperDeep  = [infraASuperBDeepR; infraADeepBSuperR];
            yValsMidDeep    = [infraAMidBDeepR;  infraADeepBMidR];
            yValsSuperSuper = infraPairCorrSuperR;
            yValsDeepDeep   = infraPairCorrDeepR;
            typeVal         = 'Infraslow';
    end

    figure;%('position',[932,103,709,838]);
    % idx = 1; 
    % for iBand = 1:5 % Feedforward
    %     subplot(2,5,idx); showExpFit(distValsNew,yValsSuperMid(:,iBand),12,-0.5,-0.4);
    %     title([bandLabels{iBand} ' - Super x Middle ']); xlim([0 20]); ylim([-0.6 1]); axis square;
    % 
    %     subplot(2,5,idx+5); showExpFit(distValsNew,yValsMidDeep(:,iBand),12,-0.5,-0.4);
    %     title([bandLabels{iBand} ' - Deep  x Middle ']);xlim([0 20]); ylim([-0.6 1]);axis square;
    % 
    %     if iBand==1; xlabel('Functional connectivity'); ylabel('Pairwise correlations');end
    %     sgtitle(['Feedforward: ' typeVal]);
    %     % [rhoValSM(iPlot,iBand),pValSM(iType,iBand)] = partialcorr(connValsNew,yValsSuperMid(:,iBand),distValsR);
    %     % [rhoValSD(iPlot,iBand),pValSD(iType,iBand)] = partialcorr(connValsNew,yValsSuperDeep(:,iBand),distValsR);
    %     % [rhoValMD(iPlot,iBand),pValMD(iType,iBand)] = partialcorr(connValsNew,yValsMidDeep(:,iBand),distValsR);
    %     idx = idx+1;
    %     % if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs'],'dir')
    %     %     [~,~]= mkdir(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs']);
    %     % end
    %     % f = gcf;
    %     % exportgraphics(f,['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs\' bandLabels{iBand} '-' typeVal '.png'],'Resolution',300);
    %     % close gcf;
    % 
    % end

    figure; idx = 1; 
    for iBand = 1:5 % Feedback
        subplot(3,5,idx); showExpFit(distValsNew,yValsSuperDeep(:,iBand),12,-0.5,-0.4);
        title([bandLabels{iBand} ' - Super x Deep']);xlim([0 20]); ylim([-0.6 1]);axis square;

        subplot(3,5,idx+5); showExpFit(distValsR,yValsSuperSuper(:,iBand),12,-0.5,-0.4);
        title([bandLabels{iBand} ' - Super x Super ']); xlim([0 20]); ylim([-0.6 1]); axis square;

        subplot(3,5,idx+10); showExpFit(distValsR,yValsDeepDeep(:,iBand),12,-0.5,-0.4);
        title([bandLabels{iBand} ' - Deep  x Deep ']);xlim([0 20]); ylim([-0.6 1]);axis square;

        if iBand==1; xlabel('Functional connectivity'); ylabel('Pairwise correlations');end
        sgtitle(['Feedback: ' typeVal]);

        % [rhoValSM(iPlot,iBand),pValSM(iType,iBand)] = partialcorr(connValsNew,yValsSuperMid(:,iBand),distValsR);
        % [rhoValSD(iPlot,iBand),pValSD(iType,iBand)] = partialcorr(connValsNew,yValsSuperDeep(:,iBand),distValsR);
        % [rhoValMD(iPlot,iBand),pValMD(iType,iBand)] = partialcorr(connValsNew,yValsMidDeep(:,iBand),distValsR);
        idx = idx+1;
        % if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs'],'dir')
        %     [~,~]= mkdir(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs']);
        % end
        % f = gcf;
        % exportgraphics(f,['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\BetweenCompartmentFigs\' bandLabels{iBand} '-' typeVal '.png'],'Resolution',300);
        % close gcf;

    end
end

%% Multiple linear regression
%% Partial correlation calculations
for iType = 1:3
    clear var
    switch iType
        case 1
            var = medPairCorrR;
        case 2
            var = medCorrEnvelopeR;
        case 3
            var = medCorrInfraSlowR; 
    end
    for iBand = 1:5
        [rhoVal(iType,iBand),pVal(iType,iBand)] = partialcorr(connValsR,var(:,iBand),distValsR);
    end
end
%% Checking code
iDate = 6; iFile =6;
datFileNum = datFileNumAll{iDate,1};
fileNum     = datFileNum(iFile);
probe1 = single(allProbeData{fileNum,iDate}.probe1Ch);
probe2 = single(allProbeData{fileNum,iDate}.probe2Ch);

probe = probe1;

fs              = 1e3; % Sampling frequency
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;
channels = 1:size(probe,2);

[spec,timeValsSpec,freqValsSpec] = mtspecgramc(probe,[5 2],params);

figure;
imagesc(timeValsSpec,freqValsSpec,10.*log10(squeeze(spec(:,:,15)')));
set(gca,'YDir','normal'); colormap jet; clim([-20 20]);colorbar; title('Before removal')

%%
freqIdx = freqValsSpec>=65 & freqValsSpec<=85;
powTimeBin = squeeze(sum(spec(:,freqIdx,:),2));
powTimeBindB = squeeze(sum(10.*log10(abs(spec(:,freqIdx,:))),2));

% Average powers over time
powTimeBinAvg = mean(powTimeBin,1,'omitnan');

badElecThreshHigh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
badElecThreshLow  = (median(powTimeBin,2)-5*mad(powTimeBin,1,2));

badChVal = [channels(sum((powTimeBin>badElecThreshHigh),1) >= floor(0.75*size(powTimeBin,1)) |(sum((powTimeBin<badElecThreshLow),1) >= floor(0.75*size(powTimeBin,1))))];

figure; plot(timeValsSpec,powTimeBin,color = [0.5 0.5 0.5]); hold on;
plot(timeValsSpec,badElecThreshHigh,'k',linewidth =2);
plot(timeValsSpec,badElecThreshLow,'k',linewidth =2); box off;

figure; plot(timeValsSpec,powTimeBindB,color = [0.5 0.5 0.5]); hold on;
badElecThreshHighdB = (median(powTimeBindB,2)+5*mad(powTimeBindB,1,2));
badElecThreshLowdB  = (median(powTimeBindB,2)-5*mad(powTimeBindB,1,2));
plot(timeValsSpec,badElecThreshHighdB,'k',linewidth =2);
plot(timeValsSpec,badElecThreshLowdB,'k',linewidth =2); box off;
badChValdB = [channels(sum((powTimeBindB>badElecThreshHighdB),1) >= floor(0.75*size(powTimeBindB,1)) |(sum((powTimeBindB<badElecThreshLowdB),1) >= floor(0.75*size(powTimeBindB,1))))];

%% Check spectrogram before and after bad timesegment removal

iDate =5; iFile =1;
datFileNum = datFileNumAll{iDate,1};
fileNum    = datFileNum(iFile);
probe1 = single(allProbeData{fileNum,iDate}.probe1Ch);
probe2 = single(allProbeData{fileNum,iDate}.probe2Ch);

probe = probe2;
channels = 1:32;
% channels([14 22]) = [];
channels(badElecB{fileNum,iDate}) = [];
% [bG,aG] = butter(3,[65 85]./(fs/2),'bandpass');
% probeBL = envelope(single(filtfilt(bG,aG,double(probe(:,channels)))),5);

probeBL = envelope(probe(:,channels));
chLim = 15; thresh = 4; timeBin = 50;

badTimeThreshHigh = median(probeBL(:,chLim:end),'all','omitnan') +thresh*mad(probeBL(:,chLim:end),[],[1,2]);
badTimeThreshLow  = median(probeBL(:,chLim:end),'all','omitnan') -thresh*mad(probeBL(:,chLim:end),[],[1,2]);
badTimeIndOld = find(mean(probeBL(:,chLim:end),2,'omitnan')>=badTimeThreshHigh | mean(probeBL(:,chLim:end),2,'omitnan')<=badTimeThreshLow);
badTimeInd = []; badTimes = []; clear badTimesCell

badTimeInd =[(badTimeIndOld-timeBin)  (badTimeIndOld+timeBin)];
badTimesCell = arrayfun(@(x) badTimeInd(x,1):badTimeInd(x,2), 1:size(badTimeInd,1), 'UniformOutput', false);
badTimes = unique([badTimesCell{:}]);
badTimes(badTimes>size(probeBL,1)) = [];
badTimes(badTimes<=0) = [];

% timeDiff = [100 diff(badTimeIndOld')];
% loc      = find(timeDiff>=100);
% Interpolation
% pFinal = mean(probe(:,15:end),2,'omitnan');
% for iLoc = 1:size(loc,2)
%     clear bT t pB removeIdx
%     pB = mean(probe(:,15:end),2,'omitnan');
%     t = 1:size(pB,1);
%     if iLoc~=size(loc,2)
%         bT = badTimeIndOld(loc(iLoc)):badTimeIndOld(loc(iLoc+1)-1);
%     else
%         bT = badTimeIndOld(loc(iLoc)):badTimeIndOld(end);
%     end
%     if isscalar(bT)
%         removeIdx = bT-100:bT+100;
%     else
%           removeIdx = bT(1)-100:bT(end)+100;
%     end
%     t(removeIdx) = [];
%     pB(removeIdx) =[];
%     newVals = interp1(t,pB,removeIdx,'previous');
%     pFinal(removeIdx,:) = newVals;
% end
%
% figure;plot(envelope(mean(probe(:,15:end),2,'omitnan'))); hold on; plot(envelope(pFinal));
%
%  badTimes = [];
% for iLoc = 1: size(loc,2)
%     clear bT
%     if iLoc~=size(loc,2)
%         bT = badTimeIndOld(loc(iLoc)):badTimeIndOld(loc(iLoc+1)-1);
%     else
%         bT = badTimeIndOld(loc(iLoc)):badTimeIndOld(end);
%     end
%     if isscalar(bT)
%         badTimes = [badTimes (bT-100): (bT+100)];
%     else
%         badTimes = [badTimes bT(1)-100:bT(end)+100];
%     end
% end
% badTimes = unique(badTimes);
% badTimes(badTimes>size(probeBL,1)) = [];
% badTimes(badTimes<=0) = [];
pNew = probe; pNew(badTimes,:) = [];
figure;plot(envelope(mean(probe(:,chLim:end),2,'omitnan'))); hold on; plot(envelope(mean(pNew(:,chLim:end),2,'omitnan')));
yline(badTimeThreshHigh); yline(badTimeThreshLow); %plot(envelope(pNew(:,chLim:end)),color = [0.5 0.5 0.5 0.1]);

fs              = 1e3; % Sampling frequency
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;
channels = 1:size(probe,2);

[spec,timeValsSpec,freqValsSpec] = mtspecgramc(mean(probe(:,chLim:end),2,'omitnan'),[5 2],params);

[specT,timeValsSpecT,freqValsSpecT] = mtspecgramc(mean(pNew(:,chLim:end),2,'omitnan'),[5 2],params);

figure;
subplot(211);imagesc(timeValsSpec,freqValsSpec,10.*log10(squeeze(spec')));
set(gca,'YDir','normal'); colormap jet; clim([-20 20]);colorbar; title('Before removal')

subplot(212);imagesc(timeValsSpecT,freqValsSpecT,10.*log10(squeeze(specT')));
set(gca,'YDir','normal'); colormap jet; clim([-20 20]);colorbar; title('After removal');

%% Time segments method 2
clear iDate iFile dateFileNum fileNum probe1 probe2 probe spec timeValsSpec...
    freqValsSpec meanS powerMeanS badTimeThreshHigh badTimeThreshLow badTimeInd ...
    probeT badTimes
iDate = 6; iFile =3;
datFileNum = datFileNumAll{iDate,1};
fileNum    = datFileNum(iFile);
probe1 = single(allProbeData{fileNum,iDate}.probe1Ch);
probe2 = single(allProbeData{fileNum,iDate}.probe2Ch);

probe1(:,badElecA{fileNum,iDate}) = [];
probe2(:,badElecB{fileNum,iDate}) = [];

probe = probe2;


[spec,timeValsSpec,freqValsSpec] = mtspecgramc(probe,[5 2],params);
meanS = mean(10.*log10(abs(spec(:,:,15:end))),3,'omitnan');
powMeanS = squeeze(sum(meanS,2));

badTimeThreshHigh   = (median(powMeanS,1)+3.5*mad(powMeanS,1,1));
badTimeThreshLow   = (median(powMeanS,1)-3.5*mad(powMeanS,1,1));

badTimeInd = floor(timeValsSpec(powMeanS>badTimeThreshHigh | powMeanS<badTimeThreshLow)*1e3);

badTimes = [];
if ~isempty(badTimeInd)
    jumpIdx = find(isoutlier([0 diff(badTimeInd)]));
    if ~isempty(jumpIdx)
        for iL = 1:length(jumpIdx)
            startIdx = badTimeInd(jumpIdx(iL));

            if iL == length(jumpIdx)
                endIdx = badTimeInd(end);
            else
                endIdx = badTimeInd(jumpIdx(iL+1)-1);
            end

            if startIdx == endIdx
                badTimes =[badTimes startIdx-50:startIdx+50];
            else
                badTimes = [ badTimes startIdx:endIdx];
            end
        end
    else
    end
    badTimes = unique(badTimes);
end

probeT = probe; probeT(badTimes,:) = [];
% figure;plot(probe(:,20)); hold on; plot(probeT(:,20));
% yline(badTimeThreshHigh); yline(badTimeThreshLow);

figure;subplot(211); imagesc(timeValsSpec,freqValsSpec,10.*log10(squeeze(spec(:,:,15)')));
set(gca,'YDir','normal'); colormap jet; clim([-20 20]); colorbar;title('Before removal')
[specT,timeValsSpecT,freqValsSpecT] = mtspecgramc(probeT,[5 2],params);
subplot(212); imagesc(timeValsSpecT,freqValsSpecT,10.*log10(squeeze(specT(:,:,15)')));
set(gca,'YDir','normal'); colormap jet; clim([-20 20]); colorbar;title('After removal')

figure;plot(timeValsSpec,powMeanS); hold on;
yline(badTimeThreshHigh); yline(badTimeThreshLow);


%% Figures for testing

probeAGamma = filtfilt(bG,aG,probeA);
envelopeA = envelope(probeA,5);
envelopeAGamma = envelope(probeAGamma,5);

probeBGamma = filtfilt(bG,aG,probeB);
envelopeB = envelope(probeB,5);
envelopeBGamma = envelope(probeBGamma,5);

typeLabels = {'WB'; 'Gamma' ; 'Power' ; 'PowerGamma'};

for iP = 1:2
    clear probe envelopeVal envelopeGamma
    switch iP
        case 1
            probe         = probeA;
            probeGamma    = probeAGamma;
            envelopeVal   = envelopeA;
            envelopeGamma = envelopeAGamma;
        case 2
            probe         = probeB;
            probeGamma    = probeBGamma;
            envelopeVal   = envelopeB;
            envelopeGamma = envelopeBGamma;
    end

    for iPlot = 1:4
        switch iPlot
            case 1
                varPlot = probe;
            case 2
                varPlot = probeGamma;
            case 3
                varPlot = envelopeVal;
            case 4
                varPlot = envelopeGamma;
        end
        figure; subplot(121); imagesc(imgaussfilt(corr(varPlot),1));
        colormap jet; clim([0 1]); xticks(1:32); yticks(1:32); colorbar; axis image;

        subplot(122);  imagesc(mean(imgaussfilt(corr(varPlot),1),1)');
        colormap jet; clim([0 1]); xticks(1:32); yticks(1:32); colorbar;
        sgtitle(['Probe: ' num2str(iP) ' Type: ' typeLabels{iPlot} ]);

        marginalTemp = mean(imgaussfilt(corr(varPlot,'Rows','complete'),1),2,'omitnan');
        fx = abs(movmean(gradient(marginalTemp),2,'omitnan'));
        transitionCh(iP,iPlot) = find(fx == max(fx));
    end

end
%% Testing with NPMK extraction

fnameNS5 = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\04_30_2025\Electrophysiology\run08\datafile0009.ns5';

numProbe = 2; %number of probes 
chanRange = [1:32;33:64]; %channel numbers in each probes; 1 row/probe
probeName = ["probeA"; "probeB"];
freqNS5 = 30e3; 
NN32_channelMap = [18 15 17 16 22 11 21 12 31 2 29 9 32 1 20 10 30 4 19 7 27 3 25 8 28 13 23 5 26 6 24 14];

channelMap = NN32_channelMap;

OutputNS5 = openNSx (fnameNS5, 'report');

samples = OutputNS5.MetaTags.DataPoints; %recording duration (1000 samples/s)

activeChan = horzcat(channelMap, channelMap+length(channelMap));
numChan = length(activeChan);

for chanNum = 1:32
    raw1Temp(chanNum, :) = OutputNS5.Data(activeChan(chanNum),:);
end

for chanNum = 1:32
    raw2Temp(chanNum, :) = OutputNS5.Data(activeChan(chanNum+length(channelMap)),:);
end

maxDigValNS5 = OutputNS5.ElectrodesInfo.MaxDigiValue;
maxAnlgValNS5 = OutputNS5.ElectrodesInfo.MaxAnalogValue;
AnlgDigCorrectionNS5 = double(maxDigValNS5/maxAnlgValNS5);
raw1Temp(1:32,:) = raw1Temp(1:32,:)/AnlgDigCorrectionNS5;
raw2Temp(1:32,:) = raw2Temp(1:32,:)/AnlgDigCorrectionNS5;

raw1Temp = raw1Temp'; raw2Temp = raw2Temp';

raw1Temp = single(downsample(filtfilt(bH,aH,double(raw1Temp)),30));
raw2Temp = single(downsample(filtfilt(bH,aH,double(raw2Temp)),30));





%% Plot scatters for pairwise correlations for all layers
function showScatterPlotsLayers(xVal,medPairCorrSuperR,medPairCorrMidR,medPairCorrDeepR,...
    envelopePairCorrSuperR,envelopePairCorrMidR,envelopePairCorrDeepR,...
    infraPairCorrSuperR,infraPairCorrMidR,infraPairCorrDeepR,...
    xLabel,yLabel,xLim,yLim,textLocX,textLocY1,textLocY2,bandLabels)

for iBand = 1:5
    figure; plotIdx = 1;
    for iType = 1:3
        switch iType
            case 1
                plotValSuper   = medPairCorrSuperR;
                plotValMid     = medPairCorrMidR;
                plotValDeep    = medPairCorrDeepR;
                plotLabel      = 'Timeseries';
            case 2
                plotValSuper = envelopePairCorrSuperR;
                plotValMid   = envelopePairCorrMidR;
                plotValDeep  = envelopePairCorrDeepR;
                plotLabel    = 'Envelope';
            case 3
                plotValSuper   = infraPairCorrSuperR;
                plotValMid     = infraPairCorrMidR;
                plotValDeep    = infraPairCorrDeepR;
                plotLabel      = 'Infraslow';
        end


        if strcmp(xLabel,'Distance')
            subplot(3,3,plotIdx);showExpFit(xVal,plotValSuper(:,iBand),textLocX,textLocY1,textLocY2); title([plotLabel ' Superficial']); axis square;
            xlim(xLim); ylim(yLim); box off; xlabel(xLabel); ylabel(yLabel);

            subplot(3,3,plotIdx+1);showExpFit(xVal,plotValMid(:,iBand),textLocX,textLocY1,textLocY2); title([plotLabel ' Middle']); axis square;
            xlim(xLim); ylim(yLim); box off; xlabel(xLabel); ylabel(yLabel);

            subplot(3,3,plotIdx+2); showExpFit(xVal,plotValDeep(:,iBand),textLocX,textLocY1,textLocY2);title([plotLabel ' Deep']); axis square;
            xlim(xLim); ylim(yLim); box off; xlabel(xLabel); ylabel(yLabel);
        else
           subplot(3,3,plotIdx);showLinearFit(xVal,plotValSuper(:,iBand),textLocX,textLocY1,textLocY2); title([plotLabel ' Superficial']); axis square;
            xlim(xLim); ylim(yLim); box off; xlabel(xLabel); ylabel(yLabel);

            subplot(3,3,plotIdx+1);showLinearFit(xVal,plotValMid(:,iBand),textLocX,textLocY1,textLocY2); title([plotLabel ' Middle']); axis square;
            xlim(xLim); ylim(yLim); box off; xlabel(xLabel); ylabel(yLabel);

            subplot(3,3,plotIdx+2); showLinearFit(xVal,plotValDeep(:,iBand),textLocX,textLocY1,textLocY2);title([plotLabel ' Deep']); axis square;
            xlim(xLim); ylim(yLim); box off; xlabel(xLabel); ylabel(yLabel);
        end
        plotIdx = plotIdx+3;
        sgtitle([plotLabel '-' bandLabels{iBand}]);
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