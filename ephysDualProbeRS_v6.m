%% ephysDualProbeRS_v5
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
        goodRuns = [1 1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN; ...
            NaN NaN NaN NaN 1 1 1 1 1 1 NaN NaN NaN NaN NaN NaN; ...
            NaN NaN NaN NaN 1 1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN; ...
            1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN NaN; ...
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

%% Get the transition channels for the recordings...
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
        if strcmp(expDate,'09_19_2022') && iRun == 4;  continue; end
        if strcmp(expDate,'02_07_2023') && iRun == 2;  continue; end
        if strcmp(expDate,'04_11_2023') && iRun == 10; continue; end
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

%% Separate LFP into different frequencies and correlate between the electrodes
clear meanPSDEEG medPSDEEG  specAMeanAll specBMeanAll specAMedAll specBMedAll meanCorrInfraSlow medCorrInfraSlow
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma';'Spiking'};

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

        if chA(1) == 0 || chA(2)== 0 || isempty(probeA) || isempty(probeB)
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

            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||...
                    (strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))...
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

            % Get pairwise correlations between the two probes...
            maxPairCorr(iDate,iRun,iBand)  = max(corr(xA,xB),[],'all','omitnan');
            meanPairCorr(iDate,iRun,iBand) = mean(corr(xA,xB),'all','omitnan');
            medPairCorr(iDate,iRun,iBand)  = corr(median(xA,2,'omitnan'),median(xB,2,'omitnan'));%median(corr(xA,xB),'all','omitnan');%

            % Get instantaneous power and correlate the powers
            envelopeABandLimited = envelope(abs(xA),5);
            envelopeBBandLimited = envelope(abs(xB),5);

            % Within probe correlations - Envelope
            envelopeIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(envelopeABandLimited,envelopeABandLimited))),'all','omitnan');
            envelopeIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(envelopeBBandLimited,envelopeBBandLimited))),'all','omitnan');

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

            % Correlate the time series
            medPairCorrSuper(iDate,iRun,iBand) = median(corr(timeASuper,timeBSuper,'rows','complete'),'all','omitnan');
            medPairCorrMid(iDate,iRun,iBand)   = median(corr(timeAMid,timeBMid,'rows','complete'),'all','omitnan'); 
            medPairCorrDeep(iDate,iRun,iBand)  = median(corr(timeADeep,timeBDeep,'rows','complete'),'all','omitnan');
           
            envelopePairCorrSuper(iDate,iRun,iBand) = median(corr(envelopeASuper,envelopeBSuper,'rows','complete'),'all','omitnan');
            envelopePairCorrMid(iDate,iRun,iBand)   = median(corr(envelopeAMid,envelopeBMid,'rows','complete'),'all','omitnan');
            envelopePairCorrDeep(iDate,iRun,iBand)  = median(corr(envelopeADeep,envelopeBDeep,'rows','complete'),'all','omitnan');
            
            infraPairCorrSuper(iDate,iRun,iBand) = median(corr(infraASuper,infraBSuper,'rows','complete'),'all','omitnan');
            infraPairCorrMid(iDate,iRun,iBand)   = median(corr(infraAMid,infraBMid,'rows','complete'),'all','omitnan'); 
            infraPairCorrDeep(iDate,iRun,iBand)  = median(corr(infraADeep,infraBDeep,'rows','complete'),'all','omitnan');

        end
    end
end


%% Plot FC versus pairwise correlations
matSize   = size(connValsAll);
connValsR = reshape(connValsAll',[matSize(1)*matSize(2) 1]);
distValsR = reshape(distSitesAll',[matSize(1)*matSize(2) 1]);

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


nanVals       = (isnan(connValsR) | isnan(medPairCorrR(:,1)));
lowCorrVals   = (medIntraCorrAR(:,4)<=0.2 | medIntraCorrBR(:,4)<=0.2);
% oneVals       = abs(infraIntraCorrBR(:,1)-1)<0.005;
removeDataIdx = nanVals | lowCorrVals;% | oneVals;

connValsR(removeDataIdx)       = [];
distValsR(removeDataIdx)       = [];
medPairCorrR(removeDataIdx,:)  = [];
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

%% Plotting
% Pairwise correlations vs functional connectivity
showScatterPlots(connValsR,medPairCorrR,medCorrEnvelopeR,medCorrInfraSlowR,...
    'Functional Connectivity','Pairwise correlations',[-0.6 1],[-0.6 1],...
    0.8,-0.5,-0.4,bandLabels)

% Pairwise correlations vs distance
showScatterPlots(distValsR,medPairCorrR,medCorrEnvelopeR,medCorrInfraSlowR,'Distance',...
    'Pairwise correlations',[0 20],[-0.6 1],12,-0.5,-0.4,bandLabels)

% Distributions of pairwise correlations and distance
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
            ylim([0 150]); yticks(0:20:150);
        % elseif 
        %     ylim([0 100]); yticks(0:10:100);
        else
             ylim([0 60]); yticks(0:10:60);
        end      
        title([plotLabel ' - ' bandLabels{iBand}]); box off;
        plotIdxConn = plotIdxConn+1;
    end
end
%% Plotting superficial, middle, deep correlations

showScatterPlotsLayers(connValsR,medPairCorrSuperR,medPairCorrMidR,medPairCorrDeepR,...
    envelopePairCorrSuperR,envelopePairCorrMidR,envelopePairCorrDeepR,...
    infraPairCorrSuperR,infraPairCorrMidR,infraPairCorrDeepR,...
    'Functional Connectivity','Pairwise correlations',[-0.6 1],[-0.6 1],...
    0.8,-0.5,-0.4,bandLabels)
%%
showScatterPlotsLayers(distValsR,medPairCorrSuperR,medPairCorrMidR,medPairCorrDeepR,...
    envelopePairCorrSuperR,envelopePairCorrMidR,envelopePairCorrDeepR,...
    infraPairCorrSuperR,infraPairCorrMidR,infraPairCorrDeepR,...
    'Distance','Pairwise correlations',[0 20],[-0.6 1],12,-0.5,-0.4,bandLabels)

%%

for iType = 1:3
    switch iType
        case 1
            plotValSuper   = medPairCorrSuperR;
            plotValMid     = medPairCorrMidR;
            plotValDeep    = medPairCorrDeepR;
            plotLabel = 'Timeseries';
        case 2
            plotValSuper = envelopePairCorrSuperR;
            plotValMid   = envelopePairCorrMidR;
            plotValDeep  = envelopePairCorrDeepR;
            plotLabel = 'Envelope';
        case 3
            plotValSuper   = infraPairCorrSuperR;
            plotValMid     = infraPairCorrMidR;
            plotValDeep    = infraPairCorrDeepR;
            plotLabel = 'Infraslow';
    end

    for iBand = 4%1:5
        figure;
        subplot(3,1,1);
        showLinearFit(connValsR,plotValSuper(:,iBand),0.8,-0.5,-0.4); title('Superficial'); axis square;
        xlim([-0.6 1]); ylim([-0.6 1]); box off; xlabel('Functional connectivity'); ylabel('Pairwise correlations');


        subplot(312);showLinearFit(connValsR,plotValMid(:,iBand),0.8,-0.5,-0.4); title('Middle'); axis square;
        xlim([-0.6 1]); ylim([-0.6 1]); box off; xlabel('Functional connectivity'); ylabel('Pairwise correlations');


        subplot(313); showLinearFit(connValsR,plotValDeep(:,iBand),0.8,-0.5,-0.4);title('Deep'); axis square;
        xlim([-0.6 1]); ylim([-0.6 1]); box off;
        xlabel('Functional connectivity'); ylabel('Pairwise correlations');
        sgtitle([plotLabel '-' bandLabels{iBand}]);       
    end


end
% sgtitle(figTitle);


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
modelfun = @(b,x) b(1) * exp(-b(2).*x);
beta0    = [10 2];

mdl = fitnlm(xVal,yVal, modelfun, beta0);
X   = linspace(min(xVal),max(xVal),1000);
coefficients = mdl.Coefficients{:, 'Estimate'};

yFitted = coefficients(1) * exp(-coefficients(2).*X) ;
plot(X,yFitted, '-k', 'LineWidth',1);
text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end

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
        end
        xlim(xLim); ylim(yLim); box off;
        xlabel(xLabel); ylabel(yLabel);
        title([plotLabel '-' bandLabels{iBand}]);
        plotIdx = plotIdx+1;
    end
end

end

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