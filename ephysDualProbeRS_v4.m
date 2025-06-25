%% ephysDualProbeRS_v4
% This program performs analysis on simultaneously recorded LFP data from two linear electrode arrays. 
% June 19,2023 - Keerthana Manikandan
% Refer ephysDualProbeRS_v1,v2,v3 for previous iterations of analyses performed. 
% This program performs the following -
% 1. Initialize all relevant variables and parameters
% 2. Store/Retrieve the LFP (1-250 Hz) and also remove 1-5 Hz
% 3. Obtain functional connectivity, distance, heart rates, anesthesia values 
% 4. Determine and remove bad time segments and channels from the two probes 
% 5. Determine the transition channel 
% 6. Determine within probe, pairwise correlations, spectrograms and LFP powers
% 7. Plot median pairwise correlations with FC for different bands
% 8. Plot median pairwise correlations with distance for different bands
% 9. Benchmarking the fit (for FC vs Pairwise) by shuffling Probe A/ Probe B time series
% 10. Plot Average LFP Power vs Anesthesia, LFP Power vs Heart Rate, LFP Power vs EEG Power for different frequencies
clear;clc;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir));
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));
clc;

%% 1. Initializing all the relevant parameters and flags
monkeyName     = 'CharlieSheen';
hemisphere     = 'Left';
saveFigureFlag = 0;
badChProbeB    = [14 22];
fs             = 1e3; % Sampling frequency
gammaBand      = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand      = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand       = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
probeBList     = [1:13 15:21 23:32];
chOutCortex    = 1:3;
chDeep         = 30:32;

disp(['Loading all required parameter names and values for ' monkeyName]); 
[allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']); 

%% 2. Store/Retrieve the LFP (1-250 Hz) and also remove 1-5 Hz
% probe1 = {}; probe2 ={}; eeg = {}; scalpEEG = {};
disp(['Storing/Retrieving LFP data for ' monkeyName]); 
for iDate =  1:size(allDates,1)
    clc; clear expdate datFileNum datFileName
    expDate = allDates(iDate,:);
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        
        % Get the name of stored file
        if strcmp(expDate,'11_01_2021') && (fileNum == 1 || fileNum == 2) % Only for Charlie Sheen
            datFileName = 'datafile_000';
        else
            datFileName = datFileNameAll{iDate,1};
        end

        if (fileNum>=10)
            datFileName = datFileName(1:end-1);
        end

        % Check if the lfp is already stored or not
        if  ~exist([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'],'file')
            [probe1{fileNum,iDate},probe2{fileNum,iDate},eeg{fileNum,iDate},scalpEEG{fileNum,iDate}] = saveLFPDualProbe(monkeyName,expDate,fileNum,datFileName,saveFolder,serverPath,fs);
            probe1 = cellfun(@single,probe1,'UniformOutput',0);
            probe2 = cellfun(@single,probe2,'UniformOutput',0);
            eeg = cellfun(@single,eeg,'UniformOutput',0);
            scalpEEG = cellfun(@single,scalpEEG,'UniformOutput',0);
        else
            % Retrieve LFP Data
            disp(['Retrieving data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            load([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
            probe1{fileNum,iDate} = probe1Ch;
            probe2{fileNum,iDate} = probe2Ch;

            % Filter the LFP to remove 1-5 Hz.
            if ~exist('fs','var'); fs = 1e3; end
            [bL,aL] = butter(3,([6 250]./(fs/2)),'bandpass');
            probe1{fileNum,iDate} = filtfilt(bL,aL,probe1{fileNum,iDate});
            probe2{fileNum,iDate} = filtfilt(bL,aL,probe2{fileNum,iDate});
            probe1 = cellfun(@single,probe1,'UniformOutput',0);
            probe2 = cellfun(@single,probe2,'UniformOutput',0);
            
            if exist('eegCh','var')
                eeg{fileNum,iDate} = single(eegCh);
            else
                eeg{fileNum,iDate} = [];
            end

            if exist('scalpEEGCh','var')
                scalpEEG{fileNum,iDate} = single(scalpEEGCh);
            else
                scalpEEG{fileNum,iDate} = [];
            end
        end
    end
end
clear probe1Ch probe2Ch eegCh 
clc; disp('Data Stored/Retrieved');

%% 3. Obtain functional connectivity, distance, heart rates, anesthesia values 
clear distSites connSites greenMapRef 
[distSites,connVals,refSites,movSites,greenMapRef] = getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);

% Get the connectivity and distance in a vector
distSitesAll = []; connValsAll = []; heartRateValsAll = []; anesthesiaValsAll = [];

for iDate =1:size(allDates,1) 
    datFileNum        = datFileNumAll{iDate,1};
    distSitesAll      = [distSitesAll;distSites{iDate,1}(datFileNum)];
    connValsAll       = [connValsAll; squeeze(mean(connVals{iDate,1}(:,:,datFileNum),[1,2],'omitnan'))];
    anesthesiaValsAll = [anesthesiaValsAll; anesthesiaLevels{iDate,1}(datFileNum)]; 
    heartRateValsAll  = [heartRateValsAll; heartRate{iDate,1}(datFileNum)];
end

disp(['Obtained/retrieved distance between probes, connectivity values, heart rate, anesthesia levels for ' monkeyName]);
%% 

if exist(['X:\Data\' monkeyName '_SqM\Left Hemisphere\' monkeyName '_vars.mat'],'file')
   load(['X:\Data\' monkeyName '_SqM\Left Hemisphere\' monkeyName '_vars.mat']);
end

%% 4. Determine and remove bad time segments and channels from the two probes 
% Determine bad time segments and bad channels
tic;
[allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);
toc;

%% Remove bad channels and times and split LFPs into different frequency bands (alpha, beta, gamma and wideband)
for iDate =1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end
        if strcmp(expDate,'04_11_2023') && iFile ==10; continue; end
        fileNum = datFileNum(iFile);

        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP

        % Remove extra channels if present
        if size(probe1{fileNum,iDate},2) == 33; probe1{fileNum,iDate}(:,1) = []; end
        if size(probe2{fileNum,iDate},2) == 33; probe2{fileNum,iDate}(:,1) = []; end

        % Remove bad channels from Probe B
        if chInCortexProbeB{iDate}(iFile)~= 1
            probe2{fileNum,iDate}(:,badChProbeB)= [];
        end

        % Remove bad channels from both probes
        if ~isempty(badElecA{fileNum,iDate})
            probe1{fileNum,iDate}(:,badElecA{fileNum,iDate}) = [];
        end

        if ~isempty(badElecB{fileNum,iDate})
            probe2{fileNum,iDate}(:,ismember(probeBList,badElecB{fileNum,iDate})) = [];
        end

        % Remove the bad time segments
        if ~isempty(allBadTimes{fileNum,iDate})
            probe1{fileNum,iDate}(allBadTimes{fileNum,iDate},:) = [];
            probe2{fileNum,iDate}(allBadTimes{fileNum,iDate},:) = [];
        end
    end
end 
clc; disp(['Determined bad channels, time segments for ' monkeyName]);

%% 5. Get the transition channels for the recordings...
% Eliminate bad channels and obtain the transition channels by
% computing slope of the marginals (obtained by averaging the intra-probe
% correlograms) - Mean is used here for sensitivity to outliers which is
% needed for picking the transition channel
clear estChInCortexA estChInCortexB
for iDate =1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    clc; disp(expDate);
    
    for iFile = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end
        if strcmp(expDate,'04_11_2023') && iFile ==10; continue; end 
        fileNum = datFileNum(iFile);

        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP

        % Retrieve probe based information
        clear probeA probeB marginalA marginalB fxA fxB
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate}; 

        % Get mean intra probe marginals from wideband range for Probe A 
        marginalA = mean(corr(probeA,'Rows','complete'),2,'omitnan');

        % Get the slope of the marginals for Probe A
        fxA = abs(movmean(gradient(marginalA),2,'omitnan')); 

        % Find the channel with maximum slope
        if strcmp(expDate,'11_01_2021') || strcmp(expDate,'01_01_2021')
            estChInCortexA{iDate}(fileNum,1) = 1;
        else
            estChInCortexA{iDate}(fileNum,1) = find(fxA == max(fxA)); 
        end

        if estChInCortexA{iDate}(fileNum,1)==20
            estChInCortexA{iDate}(fileNum,1) = (estChInCortexA{iDate}(fileNum,1) - 19);
        elseif estChInCortexA{iDate}(fileNum,1)>20
            estChInCortexA{iDate}(fileNum,1) = (estChInCortexA{iDate}(fileNum,1) - 20);
        end

        estChInCortexA{iDate}(fileNum,2) = size(probeA,2);

        % Repeat the above process for Probe B
        % Get mean marginals 
        marginalB = mean(corr(probeB,'Rows','complete'),2,'omitnan');
        fxB = abs(movmean(gradient(marginalB),2,'omitnan')); % Get the slope of marginals

        % Find the channel with maximum slope
        estChInCortexB{iDate}(fileNum,1) = find(fxB == max(fxB));

        if estChInCortexB{iDate}(fileNum,1) == 20
            estChInCortexB{iDate}(fileNum,1) = (estChInCortexB{iDate}(fileNum,1) - 19);
        elseif estChInCortexB{iDate}(fileNum,1)>20
            estChInCortexB{iDate}(fileNum,1) = (estChInCortexB{iDate}(fileNum,1) - 20);
        end
        estChInCortexB{iDate}(fileNum,2) = size(probeB,2);
    end
end 
clc; disp(['Determined transition channels for ' monkeyName]);

%% 6. Determine within probe, pairwise correlations, spectrograms and LFP powers
rowIdx = 1; eegFlag = []; 
clear meanPSDEEG medPSDEEG  specAMeanAll specBMeanAll specAMedAll specBMedAll meanCorrInfraSlow medCorrInfraSlow
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

[z,p,k] = cheby1(2,1,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

bandLabels = {'Alpha band'; 'Beta band'; 'Gamma band';'Wideband'}; 
tic;

for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];

    %     if ~exist([saveFolder '\pairCorrVars.mat'],'file')
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        clear probeA probeB

        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if ~isempty(eeg{fileNum,iDate}); eegFlag(rowIdx) = 1; else eegFlag(rowIdx) = 0; end

        if isempty(probeA) || isempty(probeB)
            meanPairCorr(rowIdx,:) = NaN; maxPairCorr(rowIdx,:)  = NaN; medPairCorr(rowIdx,:)  = NaN;

            meanIntraCorrA(rowIdx,:) = NaN; meanIntraCorrB(rowIdx,:) = NaN;
            medIntraCorrA(rowIdx,:)  = NaN; medIntraCorrB(rowIdx,:)  = NaN;

            rowIdx = rowIdx+1;
            continue;
        end

        % Get pairwise correlations
        corrA = max(imgaussfilt(corr(probeA),1),0);
        corrB = max(imgaussfilt(corr(probeB),1),0);
        lowIntraCorr(rowIdx,1) = mean(corrA,'all')<0.29|  mean(corrB,'all')<0.29;

        % Get mean, median and maximum pairwise correlations for different
        % frequency bands...
        clear timeSeriesCorr sizeSpecA sizeSpecB specA_R specB_R freqSpecCorr
        for iBand = 1:4
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||...
                    (strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))...
                    ||(strcmp(expDate,'09_19_2022') && (fileNum == 4)) || (isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate}))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||

                meanPairCorr(rowIdx,iBand)   = NaN; maxPairCorr(rowIdx,iBand)    = NaN; medPairCorr(rowIdx,iBand)  = NaN;
                meanIntraCorrA(rowIdx,iBand) = NaN; meanIntraCorrB(rowIdx,iBand) = NaN;
                medIntraCorrA(rowIdx,iBand)  = NaN; medIntraCorrB(rowIdx,iBand)  = NaN;
                continue;
            end

            clear xA yA chA chB
            chA(1) = estChInCortexA{iDate}(fileNum,1);

            if chA(1) == 0
                meanPairCorr(rowIdx,iBand)   = NaN; maxPairCorr(rowIdx,iBand)    = NaN;
                medPairCorr(rowIdx,iBand)    = NaN; meanIntraCorrA(rowIdx,iBand) = NaN;
                meanIntraCorrB(rowIdx,iBand) = NaN; medIntraCorrA(rowIdx,iBand)  = NaN;
                medIntraCorrB(rowIdx,iBand)  = NaN;
                continue;
            end
            chA(2) = size(probeA,2);

            chB(1) = estChInCortexB{iDate}(fileNum,1);

            if chB(1) == 0
                meanPairCorr(rowIdx,iBand) = NaN;maxPairCorr(rowIdx,iBand) = NaN;
                medPairCorr(rowIdx,iBand) = NaN;
                continue;
            end

            if ~strcmp(expDate,'08_08_2022')
                chB(2) = size(probeB,2);
            else
                chB(2) = chB(1);
            end

            % Get spectrograms for LFP (channels inside cortex) and EEG
            clear specA specB specEEG timeValsSpec freqValsSpec eegGood envelopeABandLimited envelopeBBandLimited infraSlowA infraSlowB
            [specA,timeValsSpec,freqValsSpec] = mtspecgramc(probeA,[5 2],params);
            [specB,~,~]                       = mtspecgramc(probeB,[5 2],params);

            if ~isempty(eeg{fileNum,iDate})
                eegGood                             = double(eeg{fileNum,iDate});
                eegGood(allBadTimes{fileNum,iDate}) = [];
                [eegPSD,psdFreq]                    = mtspectrumc(eegGood,params);

            end
            probeA = double(probeA);
            probeB = double(probeB);
            switch iBand % Get the correlation vs connectivity vs distance plots for different bands...

                case 1 % Alpha band
                    xA      = filtfilt(bA,aA,probeA(:,chA(1):chA(2)));
                    yA      = filtfilt(bA,aA,probeB(:,chB(1):chB(2)));
                    fInd    = freqValsSpec>=8 & freqValsSpec<=12; % Alpha
                    if ~isempty(eeg{fileNum,iDate})
                        fIndEEG = psdFreq>=8 & psdFreq<=12;
                    end

                case 2 % Beta band
                    xA      = filtfilt(bB,aB,probeA(:,chA(1):chA(2)));
                    yA      = filtfilt(bB,aB,probeB(:,chB(1):chB(2)));
                    fInd    = freqValsSpec>=13 & freqValsSpec<=30; % Beta
                    if ~isempty(eeg{fileNum,iDate})
                        fIndEEG = psdFreq>=13 & psdFreq<=30;
                    end

                case 3% Gamma band
                    xA      = filtfilt(bG,aG,probeA(:,chA(1):chA(2)));
                    yA      = filtfilt(bG,aG,probeB(:,chB(1):chB(2)));
                    fInd    = freqValsSpec>=30 & freqValsSpec<=90; % Gamma
                    if ~isempty(eeg{fileNum,iDate})
                        fIndEEG = psdFreq>=30 & psdFreq<= 90;
                    end

                case 4 % Wideband
                    xA      = probeA(:,chA(1):chA(2));
                    yA      = probeB(:,chB(1):chB(2));
                    fInd    = true(1,length(freqValsSpec)); % Wideband
                    if ~isempty(eeg{fileNum,iDate})
                        fIndEEG = true(1,length(psdFreq));
                    end
            end
              numChA = size(xA,2); 
              numChB = size(yA,2);
              
            % Get the intra probe correlations for channels inside the
            % cortex...
            meanIntraCorrA(rowIdx,iBand) = mean(corr(xA,'rows','complete'),'all','omitnan');
            meanIntraCorrB(rowIdx,iBand) = mean(corr(yA,'rows','complete'),'all','omitnan');
            medIntraCorrA(rowIdx,iBand)  = median(corr(xA,'rows','complete'),'all','omitnan');
            medIntraCorrB(rowIdx,iBand)  = median(corr(yA,'rows','complete'),'all','omitnan');

            % Get pairwise correlations between the two probes...
            maxPairCorr(rowIdx,iBand)  = max(corr(xA,yA),[],'all','omitnan');
            meanPairCorr(rowIdx,iBand) = mean(corr(xA,yA),'all','omitnan');
            medPairCorr(rowIdx,iBand)  = median(corr(xA,yA),'all','omitnan');

            if isempty(timeValsSpec); continue; end

            %             Get spectrogram and LFP powers for different probes
            %             Get the mean/median powers across time and frequencies -
            %             All channels....
            meanSpecA{fileNum,iDate}(:,iBand) = squeeze(mean(10.*log10(abs(specA(:,fInd,:))),[1 2],'omitnan'));
            meanSpecB{fileNum,iDate}(:,iBand) = squeeze(mean(10.*log10(abs(specB(:,fInd,:))),[1 2],'omitnan'));
            medSpecA{fileNum,iDate}(:,iBand)  = squeeze(median(10.*log10(abs(specA(:,fInd,:))),[1 2],'omitnan'));
            medSpecB{fileNum,iDate}(:,iBand)  = squeeze(median(10.*log10(abs(specB(:,fInd,:))),[1 2],'omitnan'));

            % Get overall powers from spectrogram for channels inside
            % cortex...
            specAMeanAll(rowIdx,iBand) = mean(meanSpecA{fileNum,iDate}(chA(1):chA(2),iBand),'omitnan');
            specBMeanAll(rowIdx,iBand) = mean(meanSpecB{fileNum,iDate}(chB(1):chB(2),iBand),'omitnan');
            specAMedAll(rowIdx,iBand)  = mean(medSpecA{fileNum,iDate}(chA(1):chA(2),iBand),'omitnan');
            specBMedAll(rowIdx,iBand)  = mean(medSpecB{fileNum,iDate}(chB(1):chB(2),iBand),'omitnan');

%             % Get the infraslow LFP 
%             infraSlowTimeA = filtfilt(sos,g, double([xA(1:1e3,:); xA; xA(1:1e3,:)])); 
%             infraSlowTimeB = filtfilt(sos,g, double([yA(1:1e3,:); yA; yA(1:1e3,:)])); 
%             infraSlowTimeA = single(infraSlowTimeA(1e3+1:(end-1e3),:));
%             infraSlowTimeB = single(infraSlowTimeB(1e3+1:(end-1e3),:)); 
% 
%             meanCorrInfraSlowLFP(rowIdx,iBand) = mean(corr(infraSlowTimeA,infraSlowTimeB),'all','omitnan'); 
%             medCorrInfraSlowLFP(rowIdx,iBand)  = median(corr(infraSlowTimeA, infraSlowTimeB),'all','omitnan'); 

            % Get instantaneous power and correlate the powers
            probeABandLimited = abs(xA); % Rectifying the signal
            probeBBandLimited = abs(yA);

            % Get the envelope of the signal
            envelopeABandLimited = envelope(probeABandLimited);
            envelopeBBandLimited = envelope(probeBBandLimited);

            meanCorrEnvelope(rowIdx,iBand) =  mean(corr(envelopeABandLimited,envelopeBBandLimited),'all','omitnan');
            medCorrEnvelope(rowIdx,iBand)  =  median(corr(envelopeABandLimited,envelopeBBandLimited),'all','omitnan');

            % Filter envelope from 0.01-0.1 Hz
            envelopeABandLimited = [envelopeABandLimited(1:1e3,:); envelopeABandLimited ;envelopeABandLimited(1:1e3,:) ];
            infraSlowA = filtfilt(sos,g,double(envelopeABandLimited));
            infraSlowA = single(infraSlowA(1e3+1:(end-1e3),:));

            envelopeBBandLimited = [envelopeBBandLimited(1:1e3,:); envelopeBBandLimited ;envelopeBBandLimited(1:1e3,:)];
            infraSlowB = filtfilt(sos,g,double(envelopeBBandLimited));
            infraSlowB = single(infraSlowB(1e3+1:(end-1e3),:));

            % Correlate infraslow flucutuations in instantaneous band power
            meanCorrInfraSlow(rowIdx,iBand) =  mean(corr(infraSlowA,infraSlowB),'all','omitnan');
            medCorrInfraSlow(rowIdx,iBand)  =  median(corr(infraSlowA,infraSlowB),'all','omitnan');

            % Correlate superficial and deep channels separately... 
            
            if mod(numChA,2) == 0; sChA = floor(numChA/2);else sChA = floor((numChA+1)/2); end
            if mod(numChB,2) == 0; sChB = floor(numChB/2);else sChB = floor((numChB+1)/2); end

            superMeanCorr(rowIdx,iBand) = mean(corr(xA(:,1:sChA),yA(:,1:sChB)),'all','omitnan'); 
            superMedCorr(rowIdx,iBand)  = median(corr(xA(:,1:sChA),yA(:,1:sChB)),'all','omitnan');

            superMeanEnvelopeCorr(rowIdx,iBand) = mean(corr(envelopeABandLimited(:,1:sChA),envelopeBBandLimited(:,1:sChB)),'all','omitnan');
            superMedEnvelopeCorr(rowIdx,iBand)  = median(corr(envelopeABandLimited(:,1:sChA),envelopeBBandLimited(:,1:sChB)),'all','omitnan');

            superMeanInfraSlowCorr(rowIdx,iBand) = mean(corr(infraSlowA(:,1:sChA),infraSlowB(:,1:sChB)),'all','omitnan');
            superMedInfraSlowCorr(rowIdx,iBand)  = median(corr(infraSlowA(:,1:sChA),infraSlowB(:,1:sChB)),'all','omitnan');

            if ~isempty(median(corr(xA(:,sChA+1:end),yA(:,sChB+1:end)),'all','omitnan'))
                deepMedCorr(rowIdx,iBand)  = median(corr(xA(:,sChA+1:end),yA(:,sChB+1:end)),'all','omitnan');
                deepMeanCorr(rowIdx,iBand) = mean(corr(xA(:,sChA+1:end),yA(:,sChB+1:end)),'all','omitnan');
                deepMeanEnvelopeCorr(rowIdx,iBand) = mean(corr(envelopeABandLimited(:,sChA+1:end),envelopeBBandLimited(:,sChB+1:end)),'all','omitnan');
                deepMedEnvelopeCorr(rowIdx,iBand)  = median(corr(envelopeABandLimited(:,sChA+1:end),envelopeBBandLimited(:,sChB+1:end)),'all','omitnan');
                deepMeanInfraSlowCorr(rowIdx,iBand) = mean(corr(infraSlowA(:,sChA+1:end),infraSlowB(:,sChB+1:end)),'all','omitnan');
                deepMedInfraSlowCorr(rowIdx,iBand)  = median(corr(infraSlowA(:,sChA+1:end),infraSlowB(:,sChB+1:end)),'all','omitnan');

            else
                deepMedCorr(rowIdx,iBand)  = superMeanCorr(rowIdx,iBand);
                deepMeanCorr(rowIdx,iBand) = superMedCorr(rowIdx,iBand);
                deepMeanEnvelopeCorr(rowIdx,iBand) = superMeanEnvelopeCorr(rowIdx,iBand);
                deepMedEnvelopeCorr(rowIdx,iBand)  = superMedEnvelopeCorr(rowIdx,iBand);
                deepMeanInfraSlowCorr(rowIdx,iBand) = superMeanInfraSlowCorr(rowIdx,iBand);
                deepMedInfraSlowCorr(rowIdx,iBand)  = superMedInfraSlowCorr(rowIdx,iBand);
            end

            if exist('eegGood','var')
                meanPSDEEG(rowIdx,iBand) = squeeze(mean(10.*log10(abs(eegPSD(fIndEEG))),'omitnan'));
                medPSDEEG(rowIdx,iBand)  = squeeze(median(10.*log10(abs(eegPSD(fIndEEG))),'omitnan'));
            end
        end
        rowIdx = rowIdx+1;
    end
end
% end
toc;
clear probeA probeB xA yA  probeBBandLimited probeABandLimited
% Remove bad data points from matrices...
nanIdx = find(isnan(meanPairCorr(:,1)));

connValsAll(nanIdx)=[];
distSitesAll(nanIdx) = [];
heartRateValsAll(nanIdx) = [];
anesthesiaValsAll(nanIdx) = [];
meanIntraCorrA(nanIdx,:)=[];
meanIntraCorrB(nanIdx,:)=[];
medIntraCorrA(nanIdx,:)=[];
medIntraCorrB(nanIdx,:)=[];

meanPairCorr(nanIdx,:) = [];
medPairCorr(nanIdx,:)  = [];
maxPairCorr(nanIdx,:)  = [];

specAMeanAll(nanIdx,:) = []; specBMeanAll(nanIdx,:) = [];
specAMedAll(nanIdx,:)  = []; specBMedAll(nanIdx,:) = [];
meanPSDEEG(nanIdx,:)  = []; medPSDEEG(nanIdx,:) = [];

meanCorrEnvelope(nanIdx,:) = []; medCorrEnvelope(nanIdx,:) = [];
meanCorrInfraSlow(nanIdx,:) = []; medCorrInfraSlow(nanIdx,:) = [];
% meanCorrInfraSlowLFP(nanIdx,:) = []; medCorrInfraSlowLFP(nanIdx,:) = []; 

superMeanCorr(nanIdx,:) = []; superMedCorr(nanIdx,:)= []; 
superMeanEnvelopeCorr(nanIdx,:) = []; superMedEnvelopeCorr(nanIdx,:)= [];
superMeanInfraSlowCorr(nanIdx,:) = []; superMedInfraSlowCorr(nanIdx,:) = []; 

 deepMedCorr(nanIdx,:)  = [];  deepMeanCorr(nanIdx,:) = [];
 deepMeanEnvelopeCorr(nanIdx,:) = [];  deepMedEnvelopeCorr(nanIdx,:)  = [];
 deepMeanInfraSlowCorr(nanIdx,:) =[]; deepMedInfraSlowCorr(nanIdx,:)  = [];

% movMeanPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:) = [];
save(['X:\Data\' monkeyName '_SqM\Left Hemisphere\' monkeyName '_vars.mat'],'allBadTimes','badElecA','badElecB','estChInCortexA','estChInCortexB',...
    'connValsAll','distSitesAll','meanIntraCorrA','meanIntraCorrB', 'medIntraCorrA', 'medIntraCorrB','meanPairCorr','medPairCorr',...
    'specAMeanAll','specBMeanAll','meanPSDEEG','medPSDEEG','meanCorrEnvelope','medCorrEnvelope','meanCorrInfraSlow','medCorrInfraSlow',...
    'meanCorrInfraSlowLFP','medCorrInfraSlowLFP','superMeanCorr','superMedCorr', 'superMeanEnvelopeCorr','superMedEnvelopeCorr',...
    'superMeanInfraSlowCorr', 'superMedInfraSlowCorr','deepMedCorr', 'deepMeanCorr','deepMeanEnvelopeCorr','deepMedEnvelopeCorr',...
    'deepMeanInfraSlowCorr','deepMedInfraSlowCorr','-append');

%% 
if strcmp(monkeyName,'CharlieSheen')
    patchInfo =[ 0 0 1 0 1 0 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1 1 1 0 1 0 1 0 1 0 1 1 1 1 0 1 0 0 0 1 1 0 1 1 0 1 1 1 0 1 1 1 0 1 1 1 ]; 
elseif strcmp(monkeyName,'Whiskey')
    patchInfo = [1 0 0 1 1 0 1 0 0 1 0 1 1 0 1 0 0 0 0 1 1 1 1 0 1 1 1 0 1 0 1 1 1 1 0 1 0 0 1 1 0 1 1 0 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 ];
end 
patchInfo(nanIdx) = []; 
save(['X:\Data\' monkeyName '_SqM\Left Hemisphere\' monkeyName '_vars.mat'],'patchInfo','-append'); 

% Classify pairs as S-S; S-M; M-M ie; sensory-sensory; sensory-motor;
% motor-motor

if strcmp(monkeyName,'CharlieSheen')
    pairClassCS = { 'SM';'SM';'SM';'MM';'MM';'SM';'SM';'SS';'SS';'SS';'SM';'SM';'SM';'MM';'MM';'SM';'MM';'SS';'SS';'SS';'SS';'SS';'SS';'SS';'MM';'MM';...
        'MM';'MM';'MM';'MM';'MM';'SS';'SS';'SS';'SS';'SS';'SM';'SM';'SS';'SM';'SS';'SS';'SS';'SS';'SM';'SS';'SM';'MM';'MM';'SM';'SM';'SM';'SM';'SM';...
        'SS';'SS' };
elseif strcmp(monkeyName,'Whiskey')
    pairClassW = {'SS';'SS';'SS';'SS';'SS';'SS';'SS';'SS';'SM';'MM';'MM';'MM';'MM';'MM';'SM';'SM';'SM';'SM';'SM';'SM';'SM';'SM';'SM';'SS';'SM';'SM';...
        'SM';'SM';'SS';'SS';'SS';'SS';'SS';'SS';'SS';'MM';'MM';'MM';'MM';'MM';'MM';'SS';'SM';'SM';'SM';'SM';'SS';'SS';'SS';'SS';'SS';'MM';'SM';'MM';...
        'SM';'SM';'SM';'SS';'SS';'MM';'MM' };  
end
pairClass(nanIdx) = [];
save(['X:\Data\' monkeyName '_SqM\Left Hemisphere\' monkeyName '_vars.mat'],'pairClass','-append'); 

groups = ['SS'; 'SM'; 'MM'];

pairClass = [pairClassCS; pairClassW];
pairClass(nanIdx) = [];
smCount = find(strcmp(pairClass, 'SM'));
ssCount = find(strcmp(pairClass, 'SS'));
mmCount = find(strcmp(pairClass, 'MM')); 

%% 
iDate = 5; iFile = 3; fileNum = 3; 
[z,p,k] = cheby1(2,1,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

probeA = probe1{fileNum,iDate};
probeB = probe2{fileNum,iDate};

chA(1) = estChInCortexA{iDate}(fileNum,1);
chA(2) = size(probeA,2);
chB(1) = estChInCortexB{iDate}(fileNum,1);
chB(2) = size(probeB,2);

probeA = probeA(:,chA(1):chA(2));
probeB = probeB(:,chB(1):chB(2));

[specA,timeValsSpec,freqValsSpec] = mtspecgramc(probeA,[5 2],params);
[specB,~,~]                       = mtspecgramc(probeB,[5 2],params);

specA = 10.*log10(abs(specA));
specB = 10.*log10(abs(specB));

probeA = double(probeA); 
probeB = double(probeB); 
bandLabels = {'Alpha band'; 'Beta band'; 'Gamma band';'Wideband'}; 


%% 

figure;
plotCount =1;
for iBand = 1:3
    switch iBand 
        case 1
            pA = filtfilt(bA,aA,probeA);
            pB = filtfilt(bA,aA,probeB);
            
            fInd    = freqValsSpec>=8 & freqValsSpec<=12; % Alpha
        case 2
            pA = filtfilt(bB,aB,probeA);
            pB = filtfilt(bB,aB,probeB);
            fInd    = freqValsSpec>=13 & freqValsSpec<=30; % Beta

        case 3
            pA = filtfilt(bG,aG,probeA);
            pB = filtfilt(bG,aG,probeB);
            fInd    = freqValsSpec>=30 & freqValsSpec<=90; % Gamma
    end

    for iType = 1:3
        switch iType
            case 1
                superA  = pA(:,1:8);  superB  = pB(:,1:8);
                layer4A = pA(:,9:11); layer4B = pB(:,9:11);
                deepA   = pA(:,12:17); deepB   = pB(:,12:17);
                tLabel  = 'Time-series'; 
            case 2
                superA =  squeeze(mean(specA(:,fInd,1:8),2,'omitnan'));
                superB  = squeeze(mean(specB(:,fInd,1:8),2,'omitnan'));

                layer4A = squeeze(mean(specA(:,fInd,9:11),2,'omitnan'));
                layer4B =  squeeze(mean(specB(:,fInd,9:11),2,'omitnan'));

                deepA   = squeeze(mean(specA(:,fInd,12:17),2,'omitnan'));
                deepB   =  squeeze(mean(specB(:,fInd,12:17),2,'omitnan'));
                tLabel  = 'Power'; 
            case 3
                infraSlowA = filtfilt(sos,g,envelope(abs(pA)));
                infraSlowB = filtfilt(sos,g,envelope(abs(pB)));
                superA  = infraSlowA(:,1:8);  superB  = infraSlowB(:,1:8);
                layer4A = infraSlowA(:,9:11); layer4B = infraSlowB(:,9:11);
                deepA   = infraSlowA(:,12:17); deepB   = infraSlowB(:,12:17);
                tLabel  = 'Infraslow'; 
        end
     

    corrValsAll(1,1,iBand,iType) = max(corr(superA,superB),[],[1 2]); %mean(corr(superA,superB),[1 2]);%
    corrValsAll(1,2,iBand,iType) = max(corr(superA,layer4B),[],[1 2]);%mean(corr(superA,layer4B),[1 2]);%
    corrValsAll(1,3,iBand,iType) = max(corr(superA,deepB),[],[1 2]);%mean(corr(superA, deepB),[1 2]);%

    corrValsAll(2,1,iBand,iType) = max(corr(layer4A,superB),[],[1 2]);%mean(corr(layer4A, superB),[1 2]);%
    corrValsAll(2,2,iBand,iType) = max(corr(layer4A,layer4B),[],[1 2]); %mean(corr(layer4A, layer4B),[1 2]);%
    corrValsAll(2,3,iBand,iType) = max(corr(layer4A,deepB),[],[1 2]);%mean(corr(layer4A, deepB),[1 2]); %

    corrValsAll(3,1,iBand,iType) = max(corr(deepA,superB),[],[1 2]); %mean(corr(deepA,superB),[1 2]);%
    corrValsAll(3,2,iBand,iType) = max(corr(deepA,layer4B),[],[1 2]);%mean(corr(deepA, layer4B),[1 2]);%
    corrValsAll(3,3,iBand,iType) = max(corr(deepA,deepB),[],[1 2]);%mean(corr(deepA, deepB),[1 2]);%

  subplot(3,3,plotCount);
  imagesc(corrValsAll(:,:,iBand,iType)); colormap jet; title([tLabel ' - ' bandLabels{iBand}]); 
%   if iBand ~=3;
caxis([0 1]);% else; caxis([0 0.7]); end
  axis image; colorbar;
  plotCount = plotCount+1; 
    end


end 
 
%%
for iBand = 1:4
    switch iBand 
        case 1
            pA = filtfilt(bA,aA,probeA);
            pB = filtfilt(bA,aA,probeB);
            fInd    = freqValsSpec>=8 & freqValsSpec<=12; % Alpha
        case 2
            pA = filtfilt(bB,aB,probeA);
            pB = filtfilt(bB,aB,probeB);
            fInd    = freqValsSpec>=13 & freqValsSpec<=30; % Beta

        case 3
            pA = filtfilt(bG,aG,probeA);
            pB = filtfilt(bG,aG,probeB);
            fInd    = freqValsSpec>=30 & freqValsSpec<=90; % Gamma
        case 4
            pA = probeA;
            pB = probeB;
            fInd    = true(1,length(freqValsSpec)); % Wideband
    end 

    for iType = 1: 3
        switch iType
            case 1 
                superA  = pA(:,1:8); superACorr = corr(superA,pB(:,1:20));
                superB  = pB(:,1:8); superBCorr = corr(superB,pA(:,1:20));
                layer4A = pA(:,9:11); layer4ACorr = corr(layer4A,pB(:,1:20));
                layer4B = pB(:,9:11); layer4BCorr = corr(layer4B,pA(:,1:20));
                deepA   = pA(:,12:17); deepACorr = corr(deepA,pB(:,1:20));
                deepB   = pB(:,12:17); deepBCorr = corr(deepB,pA(:,1:20));
                tLabel = 'Time series'; 

            case 2
                superA =  squeeze(mean(specA(:,fInd,1:8),2,'omitnan'));
                superACorr = corr(superA,squeeze(mean(specB(:,fInd,1:20),2,'omitnan')));
                superB  = squeeze(mean(specB(:,fInd,1:8),2,'omitnan'));
                superBCorr = corr(superB,squeeze(mean(specA(:,fInd,1:20),2,'omitnan')));

                layer4A = squeeze(mean(specA(:,fInd,9:11),2,'omitnan'));
                layer4ACorr = corr(layer4A,squeeze(mean(specB(:,fInd,1:20),2,'omitnan')));
                layer4B =  squeeze(mean(specB(:,fInd,9:11),2,'omitnan'));
                layer4BCorr = corr(layer4B,squeeze(mean(specA(:,fInd,1:20),2,'omitnan')));

                deepA   = squeeze(mean(specA(:,fInd,12:17),2,'omitnan'));
                deepACorr = corr(deepA,squeeze(mean(specB(:,fInd,1:20),2,'omitnan')));
                deepB   =  squeeze(mean(specB(:,fInd,12:17),2,'omitnan'));
                deepBCorr = corr(deepB,squeeze(mean(specA(:,fInd,1:20),2,'omitnan')));
                tLabel = 'Power'; 

            case 3

                infraSlowA = filtfilt(sos,g,envelope(abs(pA)));
                infraSlowB = filtfilt(sos,g,envelope(abs(pB)));
                superA  = infraSlowA(:,1:8);  superB  = infraSlowB(:,1:8);
                layer4A = infraSlowA(:,9:11); layer4B = infraSlowB(:,9:11);
                deepA   = infraSlowA(:,12:17); deepB   = infraSlowB(:,12:17);

                superACorr = corr(superA,infraSlowB(:,1:20)); 
                superBCorr = corr(superB, infraSlowA(:,1:20)); 

                layer4ACorr = corr(layer4A,infraSlowB(:,1:20)); 
                layer4BCorr = corr(layer4B,infraSlowA(:,1:20)); 

                deepACorr = corr(deepA,infraSlowB(:,1:20)); 
                deepBCorr = corr(deepB,infraSlowA(:,1:20)); 
                tLabel = 'Infraslow'; 
        end 
                
 figure;
   for iP = 1:2
       subplot(1,2,iP);
       switch iP
           case 1
               superCorr = superACorr; 
               deepCorr = deepACorr;
               layer4Corr = layer4ACorr; 
               figTitle = 'Probe A Reference';
           case 2
               superCorr = superBCorr; 
               deepCorr = deepBCorr;
               layer4Corr = layer4BCorr; 
               figTitle = 'Probe B Reference';
       end
       plot(movmean(mean(superCorr,1),3),1:size(superCorr,2),'LineWidth',2); hold on;  %title(' Superficial channels'); set(gca,'YDir','Reverse'); xlim([0.1 0.8]);
%        subplot(132);
       plot(movmean(mean(layer4Corr,1),3),1:size(layer4Corr,2),'LineWidth',2);% title(' Layer 4 channels');set(gca,'YDir','Reverse'); xlim([0.1 0.8]);
%        subplot(133);
       plot(movmean(mean(deepCorr,1),3),1:size(deepCorr,2),'LineWidth',2); set(gca,'YDir','Reverse');  xlim([-0.5 1]);%title(' Deep channels');set(gca,'YDir','Reverse'); xlim([0.1 0.8]);
       title([tLabel ' - ' bandLabels{iBand} ' ' figTitle]);
       legend('Superficial','Layer 4','Deep','Location','southeast'); 
   end  
    end
end 

%% 
dateVec = [2 4 5 ];
runVec = {[1 2 3 5 ]; [1 2 3 4 5 ]; [1 2 3 4 5]};
pAlphaA = zeros(32,14);
pAlphaB = zeros(32,14);
pGammaA = zeros(32,14);
pGammaB = zeros(32,14);
col = 1; 


for iDate = 1:length(dateVec)
    for iRun = 1: length(runVec{iDate})

        probeA = probe1{runVec{iDate}(iRun),dateVec(iDate)};
        probeB = probe2{runVec{iDate}(iRun),dateVec(iDate)};

        chA(1) = estChInCortexA{dateVec(iDate)}(runVec{iDate}(iRun),1);
        chA(2) = size(probeA,2);
        chB(1) = estChInCortexB{dateVec(iDate)}(runVec{iDate}(iRun),1);
        chB(2) = size(probeB,2);

        lA = length(meanSpecA{runVec{iDate}(iRun),dateVec(iDate)}(chA(1):chA(2),1)); 
        lB = length(meanSpecB{runVec{iDate}(iRun),dateVec(iDate)}(chB(1):chB(2),1));

        pAlphaA(1:lA, col) = meanSpecA{runVec{iDate}(iRun),dateVec(iDate)}(chA(1):chA(2),1); 
        pAlphaB(1:lB, col) = meanSpecB{runVec{iDate}(iRun),dateVec(iDate)}(chB(1):chB(2),1);  

        pGammaA(1:lA, col) = meanSpecA{runVec{iDate}(iRun),dateVec(iDate)}(chA(1):chA(2),3); 
        pGammaB(1:lB, col) = meanSpecB{runVec{iDate}(iRun),dateVec(iDate)}(chB(1):chB(2),3); 
        col = col+1; 
    end
end 
pA = [pAlphaA pAlphaB];
 pG = [pGammaA pGammaB];
 xA = median(pA,2,'omitnan'); madA = mad(pA,1,2);
xG = median(pG,2,'omitnan'); madG = mad(pG,1,2);

 figure; subplot(121); plot(median(pA,2,'omitnan'),1:32); set(gca,'YDir','reverse');
 patch([(xA(1:30)-madA(1:30));  flipud((xA(1:30)+madA(1:30)))],[(1:30)' ;flipud((1:30)')],'b','FaceAlpha',0.2,'EdgeColor','none');

 subplot(122); plot(median(pG,2,'omitnan'),1:32); set(gca,'YDir','reverse');
 patch([(xG(1:30)-madG(1:30));  flipud((xG(1:30)+madG(1:30)))],[(1:30)' ;flipud((1:30)')],'b','FaceAlpha',0.2,'EdgeColor','none');


%% 7. Plot median pairwise correlations with FC
figure;
X = [ones(size(connValsAll)) connValsAll distSitesAll];
cVal = {'b';'r';'m';'k'};
for iBand = 1:4
    pairWiseCorr = double(meanCorrInfraSlow(:,iBand));
    subplot(2,2,iBand);
    scatter(connValsAll,pairWiseCorr,35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand});hold on;
    %     if ~isempty(~eegFlag)
    %         scatter(connValsAll(~eegFlag),medPairCorr(~eegFlag,iBand),35,'filled','MarkerFaceColor','r','MarkerEdgeColor',cVal{iBand});
    %     end
    
    xlim([-0.4 1]); ylim([-1 1]);
    xlabel('Functional Connectivity'); ylabel('Mean  correlations');
  hold on; box on; grid on;
    title(bandLabels{iBand});

    %     pairWiseCorr(isnan(pairWiseCorr)) = 0;
%     [b,~,~,~,statsCorr(iBand,:)] = regress(pairWiseCorr,X);
%     x1Fit = linspace(min(connValsAll),max(connValsAll),1000);
%     x2Fit = linspace(min(distSitesAll),max(distSitesAll),1000);
%     yFit = b(1) + b(2)*x1Fit + b(3)*x2Fit;

    [f,gof]  = fit(connValsAll,pairWiseCorr,'poly1');
     xFit = linspace(min(connValsAll),max(connValsAll),1000);
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    r2Final (iBand) = r2; 
    gofConn(iBand) = gof; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);
end 
%% 8. Plot median pairwise correlations with distance
figure;
for iBand = 1:4
    pairWiseCorr = double(meanPairCorr(:,iBand));

    subplot(2,2,iBand);
    scatter(distSitesAll,pairWiseCorr,35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand}); hold on;
    %      if ~isempty(~eegFlag)
    %         scatter(distSitesAll(~eegFlag),medPairCorr(~eegFlag,iBand),35,'filled','MarkerFaceColor','r','MarkerEdgeColor',cVal{iBand});
    %     end
    xlim([0 20]); ylim([-1 1]);
    xlabel('Distance (mm)'); ylabel('Mean  correlations');
    hold on; box on; grid on;
    title(bandLabels{iBand})

    xFit = linspace(min(distSitesAll),max(distSitesAll),1000);
    [f,gof,output]  = fit(distSitesAll,pairWiseCorr,'poly2');
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit.^2) +c(2)*(xFit)+c(3);
    r2Final (iBand) = r2; 
    gofDist(iBand) = gof; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);
end 

%% 10. Plot Average LFP Power vs Anesthesia, LFP Power vs Heart Rate, LFP Power vs EEG Power for different frequencies
for iPlot = 1:4
    switch iPlot
        case 1
            var = anesthesiaValsAll;
            varName = 'Anesthesia Levels (%)';
        case 2
            var = heartRateValsAll;
            varName = 'Heart rate (bpm)';

        case 3
            var = meanPSDEEG;
            varName = 'Mean EEG Power (dB)';
        case 4
            var = medPSDEEG;
            varName = 'Median EEG Power (dB)';
    end

    if iPlot >=3
        if strcmp(monkeyName,'Whiskey')
            %              var(1:4,:) = []; % Uncomment for labeling anesthesia level datapoints
            var([1:4 20:22],:) = [];
        elseif strcmp(monkeyName,'CharlieSheen')
            var(18:33,:) = []; % unreliable EEG for Charlie
        end

        zeroInd = find(var(:,4)==0);
        var(zeroInd,:)  = [];

    end
    for iLFP = 1:2
        switch iLFP
            case 1
                lfpA = specAMeanAll;
                lfpB = specBMeanAll;
                lfpTitle = 'Mean LFP Power';
            case 2
                lfpA = specAMedAll;
                lfpB = specBMedAll;
                lfpTitle = 'Median LFP Power';
        end

        figure;
        if iPlot>=3
            if strcmp(monkeyName,'Whiskey')
                %                 lfpA(1:4 ,:) = []; lfpB(1:4,:) = []; % Uncomment for labeling anesthesia level datapoints
                lfpA([1:4 20:22],:) = []; lfpB([1:4 20:22],:) = [];
            elseif strcmp(monkeyName,'CharlieSheen')
                lfpA(18:33 ,:) = []; lfpB(18:33,:) = [];
            end

            lfpA(zeroInd,:) = [];
            lfpB(zeroInd,:) = [];
        end

        for iBand = 1:4
            subplot(2,2,iBand);
            if iPlot>2
                scatter(var(:,iBand),lfpA(:,iBand),30,'red','filled'); hold on;
                scatter(var(:,iBand),lfpB(:,iBand),30,'blue','filled');
            else
                scatter(var,lfpA(:,iBand),30,'red','filled'); hold on;
                scatter(var,lfpB(:,iBand),30,'blue','filled');
            end

            % Uncomment for labeling anesthesia level datapoints
            %             if strcmp(monkeyName,'Whiskey')
            %                 scatter(var(6,iBand),lfpA(6,iBand),80,'yellow','filled'); % Burst suprression
            %                 scatter(var(7,iBand),lfpA(7,iBand),80,'green','filled');  % Deep
            %
            %                 scatter(var(6,iBand),lfpB(6,iBand),80,[0.929 0.694 0.125],'filled');
            %                 scatter(var(7,iBand),lfpB(7,iBand),80,[ 0.466 0.674 0.188],'filled');
            %             end

            title(bandLabels{iBand});
            xlabel(varName); ylabel(lfpTitle); box on; grid on;
            if iBand == 1; legend('Probe A','Probe B','Location','northeast'); end
            varTemp = var;
            aTemp = lfpA;
            bTemp = lfpB;

            if iPlot>=2
                if iPlot == 2
                    nanHR = find(isnan(varTemp));
                    varTemp(nanHR) = []; aTemp(nanHR,:) = []; bTemp(nanHR,:) = [];
                    xFit         = linspace(min(varTemp),max(varTemp),1000);
                    pairWiseCorr = aTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fA,~]       = fit(varTemp,pairWiseCorr,'poly1');
                    pairWiseCorr = bTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fB,~]       = fit(varTemp,pairWiseCorr,'poly1');
                else
                    xFit = linspace(min(varTemp(:,iBand)),max(varTemp(:,iBand)),1000);
                    pairWiseCorr = aTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fA,~]       = fit(varTemp(:,iBand),pairWiseCorr,'poly1');
                    pairWiseCorr = bTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fB,~]       = fit(varTemp(:,iBand),pairWiseCorr,'poly1');
                end

                c       = coeffvalues(fA);
                fitLine = c(1)*(xFit) +c(2);
                plot(xFit,fitLine,'Color','red','LineStyle','--','LineWidth',2);

                c       = coeffvalues(fB);
                fitLine = c(1)*(xFit) +c(2);
                plot(xFit,fitLine,'Color','blue','LineStyle','--','LineWidth',2);

            end
        end
        sgtitle([lfpTitle ' vs ' varName]);
    end
end

%% Plot within probe correlations as a function of EEG Power 
clear eegPowVals slopeEEG rmseeg envelope100 meanEEG
% get EEG power, then take specific powers
rowIdx = 1; segLen = 500; winSize = 250;
for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    clc; disp(num2str(expDate));
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        if ~isempty(eeg{fileNum,iDate})
            [powEEG{fileNum,iDate},f] = pwelch(eeg{fileNum,iDate},segLen, winSize,1:120,fs);
            for iBand = 1:6
                switch iBand
                    case 1
                        fIndEEG = f>=1 & f<4; % delta
                    case 2
                        fIndEEG = f>=8 & f<=12; % alpha
                    case 3
                        fIndEEG = f>=13 & f<=30; % beta
                    case 4
                        fIndEEG = (f>=30 & f<50) | (f>60 & f<=90); % gamma 
                    case 5
                        fIndEEG = f>=1 & f<= 30; % 1-30 Hz power
                    case 6
                        fIndEEG = f>=20 & f<=30;
                end 
                eegPowVals(rowIdx,iBand) = mean(10.*log10(abs(powEEG{fileNum,iDate}(fIndEEG))),'omitnan');
            end
            slopeEEG(rowIdx,1) = (10.*log10(abs(powEEG{fileNum,iDate}(20))) - 10.*log10(abs(powEEG{fileNum,iDate}(30))))/10;

            % Get gamma band-limited signal envelope
            gammaEEG = filtfilt(bG,aG,eeg{fileNum,iDate});
               envelopeEEG{fileNum,iDate} = movmean(envelope(gammaEEG,300),1000);
            envelope100(rowIdx,:) = envelopeEEG{fileNum,iDate}(1:3000);
            meanEEG(rowIdx,1) = mean(envelopeEEG{fileNum,iDate},'omitnan' ); % 1s moving window

            [f,p1]= getFFT(envelopeEEG{fileNum,iDate},1e3,0);
            ind = find(f-0.1956>0,1);%find(f-0.2>0,1);
            pow02Hz(rowIdx,1)= 10.*log10(p1(ind));
            % RMS of gamma bands

             rmseeg(rowIdx,1) = rms(gammaEEG);

            eegFlag(rowIdx) = 1;
        else
            eegFlag(rowIdx) = 0;
        end
        rowIdx = rowIdx+1;
    end
end
eegPowVals(nanIdx,:) = [];
slopeEEG(nanIdx)=[];
rmseeg(nanIdx)=[];
envelope100(nanIdx,:)=[];
meanEEG(nanIdx,:) =[];
pow02Hz(nanIdx)=[];
%% 
figure;
for iR = 1:16
    for iC = 1:6
        if ~isempty(envelopeEEG{iR,iC})
            plot(envelopeEEG{iR,iC}); hold on
        else
            continue;
        end
    end
end 

%%  plot all eeg
figure
for iR = 1:16
    for iC = 1:6
        if ~isempty(powEEG{iR,iC})
            plot(10.*log10(abs(powEEG{iR,iC}))); hold on
        else 
            continue;
        end
    end
end 
%% Scatter plots of EEG power and LFP intra probe correlations
eegBandLabels = {'Delta'; 'Alpha'; 'Beta'; 'Gamma'};
lfpBandLabels ={'Alpha'; 'Gamma'};
if strcmp(monkeyName,'CharlieSheen'); datSets = [1:17 34:50]; elseif strcmp(monkeyName,'Whiskey'); datSets = [1:17 23:59]; end
ind = [1 3];
for iPlot = 1:2
    figure; 
    for iS = 1:4
        subplot(2,2,iS); 
        scatter(eegPowVals(datSets,iS),meanIntraCorrA(datSets,ind(iPlot)),30,'red','filled'); hold on;
        scatter(eegPowVals(datSets,iS),meanIntraCorrB(datSets,ind(iPlot)),30,'blue','filled'); 
        xlabel(['EEG power: ' eegBandLabels{iS}]); ylabel(['LFP intraprobe correlations: ' lfpBandLabels{iPlot}]);  
        ylim([0 1]); box on; 
        xFit = linspace(min(eegPowVals(datSets,iS)),max(eegPowVals(datSets,iS)),1000);
        [fA,~]  = fit([eegPowVals(datSets,iS) ; eegPowVals(datSets,iS)] ,[meanIntraCorrA(datSets,ind(iPlot)); meanIntraCorrB(datSets,ind(iPlot))],'poly1');
        c       = coeffvalues(fA);
        fitLine = c(1)*(xFit) +c(2);
        plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);

        [fA,~]  = fit([eegPowVals(datSets,iS) ] ,[meanIntraCorrA(datSets,ind(iPlot))],'poly1');
        c       = coeffvalues(fA);
        fitLine = c(1)*(xFit) +c(2);
        plot(xFit,fitLine,'Color','r','LineStyle','--','LineWidth',2);

        [fA,~]  = fit([ eegPowVals(datSets,iS)] ,[meanIntraCorrB(datSets,ind(iPlot))],'poly1');
        c       = coeffvalues(fA);
        fitLine = c(1)*(xFit) +c(2);
        plot(xFit,fitLine,'Color','b','LineStyle','--','LineWidth',2);
    end
     sgtitle(['LFP ' lfpBandLabels{iPlot} ' intraprobe correlations vs EEG Powers']);
end 


%% The code in this block is redundant 

for iPlot =1:4
    switch iPlot
        case 1
            lfpA = meanIntraCorrA;
            lfpB = meanIntraCorrB;
            eegVals = meanPSDEEG;
            figTitle = 'Mean Intra probe correlationsvs EEG Powers';
        case 2
            lfpA = medIntraCorrA;
            lfpB = medIntraCorrB;
            eegVals = medPSDEEG;
             figTitle = 'Median Intra probe correlations vs EEG Powers';
        case 3
            pairCorr = meanPairCorr; 
            eegVals  = meanPSDEEG;
             figTitle = 'Mean pairwise vs EEG Powers';
        case 4
            pairCorr = medPairCorr; 
            eegVals = medPSDEEG; 
            figTitle = 'Median pairwise vs EEG Powers';
    end
    figure;
    if strcmp(monkeyName,'Whiskey')
         eegVals([1:4 20:22],:) = [];
        if iPlot<=2
            lfpA([1:4 20:22],:) = []; lfpB([1:4 20:22],:) = [];
        else
            pairCorr([1:4 20:22],:) = [];
        end

    elseif strcmp(monkeyName,'CharlieSheen')
        eegVals(18:33,:) = [];
        if iPlot<=2
            lfpA(18:33 ,:) = []; lfpB(18:33,:) = []; 
        else
            pairCorr(18:33,:) = [];
        end
    end
    for iBand = 1:4
        subplot(2,2,iBand);
        if iPlot <=2
            scatter(eegVals(:,iBand),lfpA(:,iBand),30,'red','filled'); hold on;
            scatter(eegVals(:,iBand),lfpB(:,iBand),30,'blue','filled');
            xlabel('EEG powers (from PSD)'); ylabel('LFP within probe correlations');
            title(bandLabels{iBand}); box on;

            xFit = linspace(min(eegVals(:,iBand)),max(eegVals(:,iBand)),1000);
            [fA,~]       = fit(eegVals(:,iBand),lfpA(:,iBand),'poly1');
            [fB,~]       = fit(eegVals(:,iBand),lfpB(:,iBand),'poly1');

            c       = coeffvalues(fA);
            fitLine = c(1)*(xFit) +c(2);
            plot(xFit,fitLine,'Color','red','LineStyle','--','LineWidth',2);

            c       = coeffvalues(fB);
            fitLine = c(1)*(xFit) +c(2);
            plot(xFit,fitLine,'Color','blue','LineStyle','--','LineWidth',2);

        else
            scatter(eegVals(:,iBand),pairCorr(:,iBand),30,'k','filled'); hold on;
            xFit = linspace(min(eegVals(:,iBand)),max(eegVals(:,iBand)),1000);
            xlabel('EEG powers (from PSD)'); ylabel('Pairwise correlations');  
            title(bandLabels{iBand}); box on; 
            [fA,~]       = fit(eegVals(:,iBand),pairCorr(:,iBand),'poly1');
            c       = coeffvalues(fA);
            fitLine = c(1)*(xFit) +c(2);
            plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);
        end
        sgtitle(figTitle);
    end
end

% 
%% Plot boxplots for each anesthesia level for both probes A and probe B
powA = []; powB = []; groupVal = []; withinA = []; withinB = []; pairCorrGroup = [];
[uniqueAnesthesia,~,index] = unique(anesthesiaValsAll);
for iL = 1:length(anesthesiaValsAll)

    if strcmp(monkeyName,'CharlieSheen')
        if anesthesiaValsAll(iL)<0.9
            groupVal = [groupVal; 1];
        else
            groupVal = [groupVal; 2];
        end
        powA = [powA; specAMedAll(iL,4)];
        powB = [powB; specBMedAll(iL,4)];
        withinA  = [withinA; meanIntraCorrA(iL,4)];
        withinB  = [withinB; meanIntraCorrB(iL,4)];
        pairCorrGroup = [pairCorrGroup; medPairCorr(iL,4)];

    elseif strcmp(monkeyName,'Whiskey') && iL>4
        if anesthesiaValsAll(iL)<=1.2
            groupVal = [groupVal; 1];
            powA     = [powA; specAMedAll(iL,4)];
            powB     = [powB; specBMedAll(iL,4)];
            withinA  = [withinA; meanIntraCorrA(iL,4)];
            withinB  = [withinB; meanIntraCorrB(iL,4)];
            pairCorrGroup = [pairCorrGroup; medPairCorr(iL,4)];

        elseif any(anesthesiaValsAll(iL)>1.2 & anesthesiaValsAll(iL)<2)
            powA     = [powA; specAMedAll(iL,4)];
            powB     = [powB; specBMedAll(iL,4)];
            withinA  = [withinA; meanIntraCorrA(iL,4)];
            withinB  = [withinB; meanIntraCorrB(iL,4)];
            pairCorrGroup = [pairCorrGroup; medPairCorr(iL,4)];
            groupVal = [groupVal; 2];
          
        end
    end
end


%% 
if strcmp(monkeyName,'CharlieSheen')
    labels = {'<0.9%','>=0.9%'};
    datVals = [1:17 34:50];
elseif strcmp(monkeyName,'Whiskey')
    labels = {'<=1.2%','>1.2%'};
    datVals = [1:15 19:length(withinA)];
end 

figure;

subplot(1,3,1); boxplot(withinA(datVals), groupVal(datVals),'Labels',labels); ylabel(' Intra probe correlations'); title( 'Probe A'); ylim([ 0 1]);
p = kruskalwallis(withinA(datVals), groupVal(datVals),'off'); text(1, 0.9,['p: ' num2str(p)]); axis square;
subplot(1,3,2); boxplot(withinB(datVals), groupVal(datVals),'Labels',labels); ylabel(' Intra probe correlations'); title( 'Probe B');ylim([ 0 1]);
p = kruskalwallis(withinB(datVals), groupVal(datVals),'off'); text(1, 0.9,['p: ' num2str(p)]);axis square;
subplot(1,3,3); boxplot([withinA(datVals);withinB(datVals)], [groupVal(datVals);groupVal(datVals)],'Labels',labels); ylabel('Intra probe correlations'); title( 'Both Probes combined');ylim([ 0 1]);
p = kruskalwallis([withinA(datVals);withinB(datVals)], [groupVal(datVals);groupVal(datVals)],'off'); text(1, 0.9,['p: ' num2str(p)]); axis square;

%%
figure; boxplot(pairCorrGroup(datVals),groupVal(datVals),'Labels',labels); ylabel(' Inter probe correlations'); ylim([ -1 1]);
p = kruskalwallis(pairCorrGroup(datVals), groupVal(datVals),'off'); text(1, 0.9,['p: ' num2str(p)]);

%% Comparison of within-probe correlations across probes A and Probes B 
figure; 
for iBand = 1:4
    subplot(2,2,iBand);
    if strcmp(monkeyName,'Whiskey')
        datPoints = [5:19 23: length(meanIntraCorrA)]; 
    elseif strcmp(monkeyName,'CharlieSheen')
        datPoints = 1:length(meanIntraCorrA);
    end 
    boxplot([meanIntraCorrA(datPoints,iBand) meanIntraCorrB(datPoints,iBand)],'Labels',{'Probe A'; 'Probe B'});
    [~,p] = ttest(meanIntraCorrA(datPoints,iBand),meanIntraCorrB(datPoints,iBand)); ylabel(' Intra probe correlations')
   ylim([0 1]); text(1,0.98,['p: ' num2str(p)]);
    title(bandLabels{iBand});
end

%%  LFP powers vs LFP correlations (pairwise and intra probe) 
figure; 
for iBand = 1:4
    subplot(2,2,iBand);
    if strcmp(monkeyName,'Whiskey')
        datPoints = [5:19 23: length(meanIntraCorrA)];
    elseif strcmp(monkeyName,'CharlieSheen')
        datPoints = 1:length(meanIntraCorrA);
    end
    scatter(meanIntraCorrA(datPoints,iBand),specAMeanAll(datPoints,iBand),30,'red','filled'); hold on;
    scatter(meanIntraCorrB(datPoints,iBand),specBMeanAll(datPoints,iBand),30,'blue','filled');
    xlim([ 0 1]); 
    if iBand == 1; ylim([ 0 20]); elseif iBand == 2; ylim([-10 10]); else; ylim([-20 0]); end 

    xFitA = linspace(min(meanIntraCorrA(datPoints,iBand)),max(meanIntraCorrA(datPoints,iBand)),1000);
    xFitB = linspace(min(meanIntraCorrB(datPoints,iBand)),max(meanIntraCorrB(datPoints,iBand)),1000);
    [fA,~]       = fit(meanIntraCorrA(datPoints,iBand),specAMeanAll(datPoints,iBand),'poly1');
    [fB,~]       = fit(meanIntraCorrB(datPoints,iBand),specBMeanAll(datPoints,iBand),'poly1');


    c       = coeffvalues(fA);
    fitLine = c(1)*(xFitA) +c(2);
    plot(xFitA,fitLine,'Color','red','LineStyle','--','LineWidth',2);

    c       = coeffvalues(fB);
    fitLine = c(1)*(xFitB) +c(2);
    plot(xFitB,fitLine,'Color','blue','LineStyle','--','LineWidth',2);
    box on; grid on;
    title(bandLabels{iBand});
    xlabel('Intra probe correlations');
    ylabel('LFP Powers');
    if iBand == 1; legend('Probe A','Probe B'); end
end
%%
figure; 
for iBand = 1:4
    subplot(2,2,iBand);
    if strcmp(monkeyName,'Whiskey')
        datPoints = [5:19 23: length(medPairCorr)];
    elseif strcmp(monkeyName,'CharlieSheen')
        datPoints = 1:length(medPairCorr);
    end
    scatter(medPairCorr(datPoints,iBand),specAMeanAll(datPoints,iBand),30,'red','filled'); hold on;
    scatter(medPairCorr(datPoints,iBand),specBMeanAll(datPoints,iBand),30,'blue','filled');
    xlim([ 0 1]); 
    if iBand == 1; ylim([ 0 20]); elseif iBand == 2; ylim([-10 10]); else; ylim([-20 0]); end 

    xFitA = linspace(min(medPairCorr(datPoints,iBand)),max(medPairCorr(datPoints,iBand)),1000);
    xFitB = linspace(min(medPairCorr(datPoints,iBand)),max(medPairCorr(datPoints,iBand)),1000);
    [fA,~]       = fit(medPairCorr(datPoints,iBand),specAMeanAll(datPoints,iBand),'poly1');
    [fB,~]       = fit(medPairCorr(datPoints,iBand),specBMeanAll(datPoints,iBand),'poly1');


    c       = coeffvalues(fA);
    fitLine = c(1)*(xFitA) +c(2);
    plot(xFitA,fitLine,'Color','red','LineStyle','--','LineWidth',2);

    c       = coeffvalues(fB);
    fitLine = c(1)*(xFitB) +c(2);
    plot(xFitB,fitLine,'Color','blue','LineStyle','--','LineWidth',2);
    box on; grid on;
    title(bandLabels{iBand});
    xlabel('Pairwise correlations');
    ylabel('LFP Powers');
    if iBand == 1; legend('Probe A','Probe B'); end
end
%% Plot the median/mean spectrograms across channels  

for iDate = 1:size(allDates,1)
    clc; clear expdate datFileNum datFileName
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        clear probeA probeB
        fileNum = datFileNum(iFile);
        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP
        for iProbe = 1:2
            switch iProbe
                case 1
                    probeValMean = meanSpecA{fileNum,iDate};
                    probeValMed  = medSpecA{fileNum,iDate};
                    transitionCh = estChInCortexA{iDate}(fileNum,1);
                    probeLabel   = 'A'; 
                case 2
                    probeValMean = meanSpecB{fileNum,iDate};
                    probeValMed  = medSpecB{fileNum,iDate};
                    transitionCh = estChInCortexB{iDate}(fileNum,1);
                    probeLabel   = 'B';
            end
            if isempty(probeValMean); continue; end 
            if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs\powProbe_' probeLabel '_file_ ' num2str(iFile) '.png'],'file') || ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs\powProbe_' probeLabel '_file_ ' num2str(iFile) '.eps'],'file')
                clc; disp(['Plotting for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);
                figure('units','normalized','outerposition',[0 0 1 1]);

                for iBand = 1:4
                    subplot(2,2,iBand);
%                     plot(probeValMean(:,iBand),1:length(probeValMean(:,iBand)),'LineWidth',1); hold on;
                    plot(movmean(probeValMed(:,iBand),2),1:length(probeValMed(:,iBand)),'LineWidth',1);
                    ylim([transitionCh length(probeValMean)]); xlabel(' Power (dB)');
                    yticks(transitionCh:length(probeValMean)); ylabel('Channels');set(gca,'YDir','reverse');
                    title([ bandLabels{iBand}]);
%                     if iBand == 1; legend('Mean Power','Median Power','Location','northeast'); end
                end
                sgtitle(strrep(['Exp Date:' expDate ' File: ' num2str(iFile) ' Probe ' probeLabel],'_','\_'));

                if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs'],'dir')
                    [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs']);
                end

                if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs\powProbe_' probeLabel '_file_ ' num2str(iFile) '.png'],'file')
                    f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs\powProbe_' probeLabel '_file_ ' num2str(iFile) '.png'],'Resolution',300);
                elseif ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs\powProbe_' probeLabel '_file_ ' num2str(iFile) '.eps'],'file')
                    f = gcf;
                    exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\LaminarLFPs\powProbe_' probeLabel '_file_ ' num2str(iFile) '.eps'],'ContentType','vector');
                end

                close gcf;
            end
        end

    end
end

%% Infra slow changes in power of bandlimited LFP 
[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

% n = fir1(30,[0.01 0.1]./(fs/2));

rowIdx  = 1; tic;
for iDate = 2:size(allDates,1)
     clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum)
        clear probeA probeB infraA infraB envelopeA envelopeB
        fileNum = datFileNum(iFile);
        clc; disp(['Obtaining slow changes in band limited LFP for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        if isempty(probeA) || isempty(probeB); rowIdx = rowIdx+1; continue; end

        % Get channels inside cortex...
        chA = estChInCortexA{iDate}(fileNum,:);
        chB = estChInCortexB{iDate}(fileNum,:);

        if chA(1)== 0 || chB(1) == 0
            rowIdx = rowIdx+1;
            continue;
        end
        
        probeA = probeA(:,chA(1):chA(2));
        probeB = probeB(:,chB(1):chB(2));

        for iBand = 1:3
            clear probeATemp probeBTemp
            switch iBand
                case 1 % Alpha band 
                     probeATemp = filtfilt(bA,aA,probeA);
                     probeBTemp = filtfilt(bA,aA,probeB);
                case 2 % Beta
                    probeATemp = filtfilt(bB,aB,probeA);
                     probeBTemp = filtfilt(bB,aB,probeB);
                case 3 % Gamma band 
                    probeATemp = filtfilt(bG,aG,probeA);
                    probeBTemp = filtfilt(bG,aG,probeB);
            end 
            
            % Rectify the signal
            probeATemp = abs(probeATemp);
            probeBTemp = abs(probeBTemp);

            % Get the envelope of the signal
            envelopeA = envelope(probeATemp);
            envelopeB = envelope(probeBTemp);

            % Filter envelope from 0.01-0.1 Hz
            infraA = filtfilt(sos,g,envelopeA); %filtfilt(n,1,envelopeA);%
            infraB = filtfilt(sos,g,envelopeB); %filtfilt(n,1,envelopeB);%

            % Get pairwise correlations
            pairCorrInfra(rowIdx,iBand) = median(corr(infraA,infraB),'all','omitnan');
        end
        rowIdx = rowIdx+1; 
    end
end 
toc;
zeroInd = find(pairCorrInfra(:,1)==0);
pairCorrInfra(zeroInd,:) = [];
% datVals = [1:15 19:length(pairCorrInfra)];


%%
datPoints = [1:19 22:length(pairCorrInfra) ];
datVals = datPoints; 
figure;
cVal = {'b';'r';'m';'k'};
for iBand = 1:3
    subplot(1,3,iBand);
    scatter(connValsAll(datPoints),pairCorrInfra(datVals,iBand),35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand});hold on; 
%     if ~isempty(~eegFlag)
%         scatter(connValsAll(~eegFlag),medPairCorr(~eegFlag,iBand),35,'filled','MarkerFaceColor','r','MarkerEdgeColor',cVal{iBand});
%     end 

    xlim([-0.4 1]); ylim([-1 1]);
    xlabel('Functional Connectivity'); ylabel('Infra slow changes in power');
    hold on; box on; grid on; axis square
    title(bandLabels{iBand})

    xFit = linspace(min(connValsAll(datPoints)),max(connValsAll(datPoints)),1000);
    pairWiseCorr = pairCorrInfra(datVals,iBand);
    pairWiseCorr(isnan(pairWiseCorr)) = 0;
    [f,gof]  = fit(connValsAll(datPoints),pairWiseCorr,'poly1');
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    r2Final (iBand) = r2; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);

end

%% Time binned pairwise correlations... time bin = 1000 ms, moving window 500 ms
rowIdx = 1;  
winLen   = 1*fs; % Bin length - 1 s
binWidth = 0.5*fs; % Sliding window length  - 0.5 s
for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if isempty(probeA) || isempty(probeB)
            corrValsBinned(rowIdx,iBand) = NaN;
            
            rowIdx = rowIdx+1;
            continue; 
        end

        if size(probeB,2) == 33; probeB(:,1) = []; end
        if size(probeA,2) == 33; probeA(:,1) = []; end

%         disp('Obtaining pairwise correlations...');
        % Get mean, median and maximum pairwise correlations for different
        % frequency bands... 
        for iBand = 1:4
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||(strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))||(strcmp(expDate,'09_19_2022') && (fileNum == 4))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||
                corrValsBinned(rowIdx,iBand) = NaN;
                continue;
            end

            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate})
                corrValsBinned(rowIdx,iBand) = NaN;
                continue;
            end

            clear xA yA chA chB
            chA = estChInCortexA{iDate}(fileNum,:);

            if chA(1) == 0
                corrValsBinned(rowIdx,iBand) = NaN;
                continue;
            end

            chB = estChInCortexB{iDate}(fileNum,:);
            if chB(1) == 0
                 corrValsBinned(rowIdx,iBand) = NaN;
                continue;
            end

            switch iBand % Get the correlation vs connectivity vs distance plots for different bands...

                case 1 % Alpha band
                    xA = filtfilt(bA,aA,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bA,aA,probeB(:,chB(1):chB(2)));

                case 2 % Beta band
                    xA = filtfilt(bB,aB,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bB,aB,probeB(:,chB(1):chB(2)));

                case 3 % Gamma band
                    xA = filtfilt(bG,aG,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bG,aG,probeB(:,chB(1):chB(2)));

                case 4 % Wideband
                    xA = probeA(:,chA(1):chA(2));
                    yA = probeB(:,chB(1):chB(2));
            end

            % Get lagged correlations for the two probes
            sizeA = size(xA,2); sizeB = size(yA,2); 
            ind = 1;
            clear corrDatAB
            L = size(xA,1);
           

            for iW = 1:binWidth:L-winLen+1
                corrDatAB(ind,:,:) = corr(xA(iW:iW+winLen-1,:),yA(iW:iW+winLen-1,:));
                ind = ind+1;
            end
         
            corrValsBinned(rowIdx,iBand) = median(corrDatAB,'all','omitnan'); % Taking the median correlations
        end 

         rowIdx = rowIdx+1;
    end
end 
%%
figure; 
for iBand = 1:4
  
    subplot(2,2,iBand);
    scatter(connValsAll,squeeze(corrValsBinned(:,iBand)),20,'filled');
    xlabel('Functional connectivity'); xlim([-0.3 1]);
    ylabel('Pairwise correlations'); ylim([-1 1]);
    hold on; box on;title(bandLabels{iBand}); grid on;

    xFit = linspace(min(connValsAll),max(connValsAll),1000);
    [f,gof]  = fit(connValsAll,squeeze(corrValsBinned(:,iBand)),'poly1');

    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);
    text(0.6,-0.8,['R^2: ' num2str(r2)]);
end 

f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\binnedCorrPlots.png'],'Resolution',300);
close gcf;

%% 9. Benchmarking the fit (for FC vs Pairwise) by shuffling Probe A/ Probe B time series
clear corrValsShuffledB corrValsShuffledA
rowIdx = 1; 
for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);

        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if isempty(probeA) || isempty(probeB); continue; end 

        if isempty(probeA) || isempty(probeB)
            corrValsShuffledA(rowIdx,iBand) = NaN;
            corrValsShuffledB(rowIdx,iBand) = NaN;
            rowIdx = rowIdx+1;
            continue; 
        end

        % Get mean, median and maximum pairwise correlations for different
        % frequency bands... 
        for iBand = 1:4
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||(strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))||...
                    (strcmp(expDate,'09_19_2022') && (fileNum == 4))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||
                corrValsShuffledA(rowIdx,iBand) = NaN;
                corrValsShuffledB(rowIdx,iBand) = NaN;
                continue;
            end

            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate})
                corrValsShuffledA(rowIdx,iBand) = NaN;
                corrValsShuffledB(rowIdx,iBand) = NaN;
                continue;
            end

            clear xA yA chA chB
            chA = estChInCortexA{iDate}(fileNum,:);

            if chA(1) == 0
                corrValsShuffledA(rowIdx,iBand) = NaN;
                corrValsShuffledB(rowIdx,iBand) = NaN;
                continue;
            end

            chB = estChInCortexB{iDate}(fileNum,:);
            if chB(1) == 0
                 corrValsShuffledA(rowIdx,iBand) = NaN;
                 corrValsShuffledB(rowIdx,iBand) = NaN;
                continue;
            end
            clear probeANew probeBNew
            for iC = 1:size(probeA,2)
                probeANew(:,iC) = shuffle(probeA(:,iC));
            end 

            for iC = 1:size(probeB,2)
                 probeBNew(:,iC) = shuffle(probeB(:,iC)); 
            end 

            switch iBand % Get the correlation vs connectivity vs distance plots for different bands...

                case 1 % Alpha band
                    xA    = filtfilt(bA,aA,probeA(:,chA(1):chA(2)));
                    yA    = filtfilt(bA,aA,probeB(:,chB(1):chB(2)));
                    xANew = filtfilt(bA,aA,probeANew(:,chA(1):chA(2)));
                    yANew = filtfilt(bA,aA,probeBNew(:,chB(1):chB(2)));

                case 2 % Beta band
                    xA    = filtfilt(bB,aB,probeA(:,chA(1):chA(2)));
                    yA    = filtfilt(bB,aB,probeB(:,chB(1):chB(2)));
                    xANew = filtfilt(bB,aB,probeANew(:,chA(1):chA(2))); 
                    yANew = filtfilt(bB,aB,probeBNew(:,chB(1):chB(2)));

                case 3 % Gamma band
                    xA    = filtfilt(bG,aG,probeA(:,chA(1):chA(2)));
                    yA    = filtfilt(bG,aG,probeB(:,chB(1):chB(2)));
                    xANew = filtfilt(bG,aG,probeANew(:,chA(1):chA(2)));
                    yANew = filtfilt(bG,aG,probeBNew(:,chB(1):chB(2)));

                case 4 % Wideband
                    xA    = probeA(:,chA(1):chA(2));
                    yA    = probeB(:,chB(1):chB(2));
                    xANew = probeANew(:,chA(1):chA(2));
                    yANew = probeBNew(:,chB(1):chB(2));
            end
            
            corrValsShuffledA(rowIdx,iBand) = mean(corr(xA,yANew),'all','omitnan'); % Taking the median shuffled correlations
            corrValsShuffledB(rowIdx,iBand) = mean(corr(xANew,yA),'all','omitnan'); % Taking the median shuffled correlations

            pairCorrTimeSeries(rowIdx,iBand) = mean(corr(xA,yA),'all','omitnan');

        end 
         rowIdx = rowIdx+1;
    end 
end
clc; disp(['Processing complete for ' monkeyName]); 

nanVals = find(isnan(corrValsShuffledA(:,1)));
corrValsShuffledA(nanVals,:)  = []; 
corrValsShuffledB(nanVals,:)  = [];
pairCorrTimeSeries(nanVals,:) = [];

if strcmp(monkeyName,'CharlieSheen')
    colorVals =[0.6350 0.0780 0.1840];
elseif strcmp(monkeyName,'Whiskey')
    colorVals = [0 0.4470 0.7410];
end

figure; 
for iBand = 1:4
    subplot(2,2,iBand);
    scatter(connValsAll,squeeze(pairCorrTimeSeries(:,iBand)),20,colorVals,'filled');
    xlabel('Functional connectivity'); xlim([-0.3 1]);
    ylabel('Pairwise correlations'); ylim([-1 1]);
    hold on; box on;title(bandLabels{iBand}); grid on;

    xFit = linspace(min(connValsAll),max(connValsAll),1000);
    [f,gof]  = fit(connValsAll,squeeze(pairCorrTimeSeries(:,iBand)),'poly1');

    c       = coeffvalues(f);
    r2      = gof.rsquare;
    fitLine = c(1)*(xFit) +c(2);

    plot(xFit,fitLine,'Color',colorVals,'LineStyle','--','LineWidth',2);
    text(0.6,-0.6-0.1*(iPr-1),['R^2: ' num2str(r2)],'Color',colorVals);
end

%% Power values dependent on sliding window size 

for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum)
        clear channelA channelB
       fileNum = datFileNum(iFile);
        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        chA    = estChInCortexA{iDate}(fileNum,:);
        chB    = estChInCortexB{iDate}(fileNum,:);
        channelA = chA(1):chA(2);
        channelB = chB(1):chB(2);
        if isempty(probeA) || isempty(probeB); continue; end  
        if (chA(1)==0 || chB(1) == 0); continue; end  

        % Get power - instantaneous, power spectrum and spectrogram
        % depending upon the segment length 
        winSize = 125; % Sliding window length
        allSeg  = 250: 250:1000; 
        for iSeg = 1:length(allSeg)
            disp(['Segment length:' num2str(allSeg(iSeg))]);
            dispstat('','init');
            for iCh = 1:length(channelA)
                [powA{iSeg,1}{fileNum,iDate}(:,:,iCh),freqVals{iSeg},timeVals{iSeg}] = spectrogram(probeA(:,channelA(iCh)),allSeg(iSeg), winSize,1:120,fs);
                
                if ~(iCh>length(channelB))
                    [powB{iSeg,1}{fileNum,iDate}(:,:,iCh),freqVals{iSeg},timeVals{iSeg}] = spectrogram(probeB(:,channelB(iCh)),allSeg(iSeg),winSize,1:120,fs);
                end
                dispstat(['Getting the power for the two probes... ' num2str((iCh/size(probeA,2))*100) '% done']);
            end
            dispstat('Getting the power for the two probes... 100% done');
        end

    end
end 
%%
% connValsAll(8) = []; %whiskey
connValsAll([2 42]) =[]; % Charlie Sheen

%% Get the average powers for different bands 

for iSeg = 1: length(allSeg)
    powSegA = powA{iSeg,1};
    powSegB = powB{iSeg,1};
    figure('units','normalized','outerposition',[0 0 1 1]);
    for iBand = 1:4
        switch iBand
            case 1
                fVals = 8:12;
                bandName = 'Alpha band';
            case 2
                fVals = 13:30;
                bandName = 'Beta band';
            case 3
                fVals = 30:90;
                bandName = 'Gamma band ';
            case 4
                fVals = 1:120;
                bandName = 'Wideband';
        end   
        clear tempA tempB ; 
        rowIdx = 1;
        for iC = 1:size(powSegA,2)
            for iR = 1:size(powSegA,1)
                if isempty(powSegA{iR,iC}); continue; end 
                tempA = mean(abs(powSegA{iR,iC}(fVals,:,:)),[1 3],'omitnan'); 
                tempB = mean(abs(powSegB{iR,iC}(fVals,:,:)),[1 3],'omitnan');
                pairCorr(iSeg,iBand,rowIdx) = corr(tempA',tempB','rows','complete');
                rowIdx = rowIdx+1;
            end 
        end 
        subplot(2,2,iBand);
        scatter(connValsAll,squeeze(pairCorr(iSeg,iBand,:)),30,'filled');
        xlabel('Functional connectivity'); xlim([-0.3 1]);
        ylabel('Pairwise correlations'); ylim([-1 1]);
        hold on; box on;title(bandName); grid on;

        xFit = linspace(min(connValsAll),max(connValsAll),1000);
        [f,gof]  = fit(connValsAll,squeeze(pairCorr(iSeg,iBand,:)),'poly1');
        c        = coeffvalues(f);
        r2       = gof.rsquare;
        fitLine   = c(1)*(xFit) +c(2);
        plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);
        text(0.8,-0.9,['R^2: ' num2str(r2)]);
    end

    sgtitle(['Segment length: ' num2str(allSeg(iSeg))]);
    f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\PowCorr_' num2str(allSeg(iSeg)) '.png'],'Resolution',300);
    close gcf;
end



%% 7. Check the state of the animal under anesthesia and remove recordings if needed 
segLen = 500; winSize = 250;
% eegWhiskey = eeg;
tic;
dispstat('','init');
dispstat('Getting power for EEGs... ');
count = 1; 
for iDate = 1:size(eeg,2)
    for iFile = 1:size(eeg,1)
        if isempty(eeg{iFile,iDate})
            count = count+1; 
            continue; 
        else
            [powEEG{iFile,iDate},f] = pwelch(eeg{iFile,iDate},segLen, winSize,1:120,fs);
            count = count+ 1;         
        end 
           dispstat(['Getting power for EEGs...' num2str((count/numel(eeg))*100) '% Complete']);
    end 
end 
dispstat('Getting power for EEGs... 100% Complete');
toc;

%%
% figure; 
for iDate = 1%:size(eeg,2)
    for iFile = 1:size(eeg,1)
        if isempty(powEEG{iFile,iDate}) || (iDate == 3) && sum(iFile == 6:8) || ((iDate == 5) && sum(iFile == 1:7))%||((iDate == 4) 
            continue; 
        else
            plot(f,10.*log10(powEEG{iFile,iDate}),'g','LineWidth',2);  hold on; 
        end 
    end
end 

xlabel('Frequency (Hz)'); ylabel('Power (dB)'); ylim([-20 50]); 


%%
iDate= 1;
for iR = 1: 11
    clear chA chB
    chA = estChInCortexA{1,iDate}(iR,:);
    chB = estChInCortexB{1,iDate}(iR,:);
    
    if all(chA == 0) || all(chB == 0)
        continue; 
    end 

    alphaPowA(iR,1)  = median(powA{iR,iDate}(8:12,chA(1):chA(2)),'all','omitnan');
    alphaPowB(iR,1)  = median(powB{iR,iDate}(8:12,chB(1):chB(2)),'all','omitnan');
    alphaEEG(iR,1)   = median(powEEG{iR,iDate}(8:12),'all','omitnan');
    
    betaPowA(iR,1)  = median(powA{iR,iDate}(13:30,chA(1):chA(2)),'all','omitnan');
    betaPowB(iR,1)  = median(powB{iR,iDate}(13:30,chB(1):chB(2)),'all','omitnan');
    betaEEG(iR,1)   = median(powEEG{iR,iDate}(13:30),'all','omitnan');

    gamma1PowA(iR,1)  = median(powA{iR,iDate}(30:55,chA(1):chA(2)),'all','omitnan');
    gamma1PowB(iR,1)  = median(powB{iR,iDate}(30:55,chB(1):chB(2)),'all','omitnan');
    gamma1EEG(iR,1)   = median(powEEG{iR,iDate}(30:55),'all','omitnan');

    gamma2PowA(iR,1)  = median(powA{iR,iDate}(65:90,chA(1):chA(2)),'all','omitnan');
    gamma2PowB(iR,1) = median(powB{iR,iDate}(65:90,chB(1):chB(2)),'all','omitnan');
    gamma2EEG(iR,1)  = median(powEEG{iR,iDate}(65:90),'all','omitnan');


    wbPowA(iR,1)  = median(powA{iR,iDate}(:,chA(1):chA(2)),'all','omitnan');
    wbPowB(iR,1)  = median(powB{iR,iDate}(:,chB(1):chB(2)),'all','omitnan');
    wbEEG(iR,1)   = median(powEEG{iR,iDate},'all','omitnan');
end

%%

for iBand = 1:5
    figure;
    clear vals titleVal 
    switch iBand
        case 1
            vals = [alphaPowA alphaPowB alphaEEG];
            titleVal = 'Alpha band';
        case 2
            vals = [betaPowA betaPowB betaEEG];
             titleVal = 'Beta band';
        case 3
            vals = [gamma1PowA gamma1PowB gamma1EEG];
             titleVal = 'Gamma1 band';
        case 4
            vals =[gamma2PowA gamma2PowB gamma2EEG];
             titleVal = 'Gamma2 band';
        case 5 
            vals = [wbPowA wbPowB wbEEG];
             titleVal = 'Wide band';
    end 
    boxplot(10.*log10(vals),'Labels',{'Probe A','Probe B','EEG'}); 
    ylabel('Power (dB)'); 
    title(titleVal);
end 

%% 


%% Shuffle a time series
 function v =shuffle(v)
     v=v(randperm(length(v)));
 end