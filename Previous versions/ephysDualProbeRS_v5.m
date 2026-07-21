%% ephysDualProbeRS_v5
% This function performs analysis on LFP recorded simultaneously from two linear electrode arrays.
% March 1, 2024 - Keerthana Manikandan
% This code performs the following for ONE MONKEY (for combined analyis,
% check this script- combinedAnalysisEphysDualProbeRS.m) - KM needs to write this scirpt now :/
% New script cos, everything got out of hand

% 1. Initialize all relevant variables and parameters
% 2. Store/Retrieve the LFP (1-250 Hz) and also remove 1-5 Hz
% 3. Obtain functional connectivity, distance, heart rates, anesthesia values
% 4. Determine and remove bad time segments and channels from the two probes
% 5. Determine the transition channel
% 6. Determine within probe, pairwise correlations, spectrograms and LFP powers
% 7. Saves all variables for one monkey

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

%% 1. Initializing all relevant variables

monkeyName     = 'Whiskey';
hemisphere     = 'Left';
saveFigureFlag = 0;
fs             = 1e3; % Sampling frequency
gammaBand      = [30 90]; [zG,pG,kG] = butter(3,gammaBand./(fs/2),'bandpass'); [sosG,gG] = zp2sos(zG,pG,kG);
alphaBand      = [8 12];  [zA,pA,kA] = butter(3,alphaBand./(fs/2),'bandpass'); [sosA,gA] = zp2sos(zA,pA,kA);
betaBand       = [13 30]; [zB,pB,kB] = butter(3,betaBand./(fs/2),'bandpass');  [sosB,gB] = zp2sos(zB,pB,kB);
% probeBList     = [1:13 15:21 23:32];
chOutCortex    = 1:3;
chDeep         = 30:32;

disp(['Loading all required parameter names and values for ' monkeyName]);
[allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate,patchInfo,pairClass] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']);

%% 2. Store/Retrieve LFP from the two probes...
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

            [probe1{fileNum,iDate},probe2{fileNum,iDate},eeg{fileNum,iDate},scalpEEG{fileNum,iDate}] = ...
                saveLFPDualProbe(monkeyName,expDate,fileNum,datFileName,saveFolder,serverPath,fs);

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
            [zL,pL,kL] = butter(3,([6 250]./(fs/2)),'bandpass');
            [sosL,gL]  = zp2sos(zL,pL,kL);

            probe1{fileNum,iDate} = filtfilt(sosL,gL,probe1{fileNum,iDate});
            probe2{fileNum,iDate} = filtfilt(sosL,gL,probe2{fileNum,iDate});

            probe1 = cellfun(@single,probe1,'UniformOutput',0);
            probe2 = cellfun(@single,probe2,'UniformOutput',0);

            if exist('eegCh','var') % Get EEG
                eeg{fileNum,iDate} = single(eegCh);
            else
                eeg{fileNum,iDate} = [];
            end

            if exist('scalpEEGCh','var') % Get Scalp EEG
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
[distSites,connVals,refSites,movSites,greenMapRef] = ...
    getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);

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

%% 4a. Determine and remove bad time segments and channels from the two probes
% Determine bad time segments and bad channels
tic;
[allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);
toc;

%% 4b. Remove bad channels and times
for iDate =1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iRun == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iRun == 2; continue; end
        if strcmp(expDate,'04_11_2023') && iRun ==10; continue; end
        fileNum = datFileNum(iRun);

        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP

        % Remove extra channels if present
        if size(probe1{fileNum,iDate},2) == 33; probe1{fileNum,iDate}(:,1) = []; end
        if size(probe2{fileNum,iDate},2) == 33; probe2{fileNum,iDate}(:,1) = []; end

        % Remove bad channels from both probes
        if ~isempty(badElecA{fileNum,iDate})
            probe1{fileNum,iDate}(:,badElecA{fileNum,iDate}) = [];
        end

        if ~isempty(badElecB{fileNum,iDate})
            probe2{fileNum,iDate}(:,badElecB{fileNum,iDate}) = [];
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

    for iRun = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iRun == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iRun == 2; continue; end
        if strcmp(expDate,'04_11_2023') && iRun ==10; continue; end
        fileNum = datFileNum(iRun);

        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP

        % Retrieve probe based information
        clear probeA probeB marginalA marginalB fxA fxB
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        % Get the transition channel as determined during the experiment
        transitionChA = chInCortexProbeA{iDate}(iRun);
        transitionChB = chInCortexProbeB{iDate}(iRun);      

        % Determine transition channel for Probe A
        if size(probeA,2)== 1
            estChInCortexA{iDate}(iRun,:) = [1 1];
        else
            % Get mean intra probe marginals from wideband range for Probe A
            marginalA = mean(corr(probeA,'Rows','complete'),2,'omitnan');

            % Get the slope of the marginals for Probe A
            fxA = abs(movmean(gradient(marginalA),2,'omitnan'));
            clear tempCh
            tempCh = find(fxA == max(fxA));

            if tempCh>=20 % Set transition to 1 if maximum slope is after channel 20
                tempCh = 1;
            end

            if abs(tempCh - transitionChA)>=10 % Set transition to the value determined from notes if difference between estimated and observed value exceeds 10
                estChInCortexA{iDate}(iRun,1) = transitionChA;

            elseif abs(tempCh - transitionChA)<= 3 % Set transition to estimated value if difference between estimated and observed value is less than or equal to 3
                estChInCortexA{iDate}(iRun,1) = tempCh;
            else
                estChInCortexA{iDate}(iRun,1) = floor((tempCh + transitionChA)/2); % Set transition to mean between estimated and observed value otherwise
            end

            if estChInCortexA{iDate}(iRun,1)~=0
                estChInCortexA{iDate}(iRun,2) = size(probeA,2);
            else
                estChInCortexA{iDate}(iRun,2) = 0; 
            end 
        end

        % Repeat the above process for Probe B      
        if size(probeB,2)== 1
            estChInCortexB{iDate}(iRun,:) = [1 1];
        else
            % Get mean marginals
            marginalB = mean(corr(probeB,'Rows','complete'),2,'omitnan');
            % Get the slope of marginals
            fxB = abs(movmean(gradient(marginalB),2,'omitnan'));
            clear tempCh
            tempCh = find(fxB == max(fxB));

            if tempCh>=20 % Set transition to 1 if max slope is after channel 20
                tempCh = 1;
            end

            if abs(tempCh - transitionChB)>=10 % Set transition to the value determined from notes if difference between estimated and observed value exceeds 10
                estChInCortexB{iDate}(iRun,1) = transitionChB;

            elseif abs(tempCh - transitionChB)<= 3 % Set transition to estimated value if difference between estimated and observed value is less than or equal to 3
                estChInCortexB{iDate}(iRun,1) = tempCh;
            else
                estChInCortexB{iDate}(iRun,1) = floor((tempCh + transitionChB)/2); % Set transition to mean between estimated and observed value otherwise
            end
            
            if estChInCortexB{iDate}(iRun,1)~=0
                estChInCortexB{iDate}(iRun,2) = size(probeB,2);
            else
                estChInCortexB{iDate}(iRun,2) = 0; 
            end

        end
    end
end
clc; disp(['Determined transition channels for ' monkeyName]);

%% 6. Calculate pairwise correlations
rowIdx = 1; eegFlag = []; 
clear meanPSDEEG medPSDEEG  specAMeanAll specBMeanAll specAMedAll specBMedAll meanCorrInfraSlow medCorrInfraSlow
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

[z,p,k] = butter(2,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

bandLabels = {'Alpha band'; 'Beta band'; 'Gamma band';'Wideband'}; 

tic;
for iDate = 1: size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        clear probeA probeB chA chB

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

        chA = estChInCortexA{iDate}(iRun,:); 
        chB = estChInCortexB{iDate}(iRun,:); 

        if chA(1) == 0 || chA(2)== 0
             rowIdx = rowIdx+1;
             meanPairCorr(rowIdx,1:4,1:2)   = NaN; maxPairCorr(rowIdx,1:4,1:2)    = NaN;
                medPairCorr(rowIdx,1:4,1:2)    = NaN; meanIntraCorrA(rowIdx,1:4,1:2) = NaN;
                meanIntraCorrB(rowIdx,1:4,1:2) = NaN; medIntraCorrA(rowIdx,1:4,1:2)  = NaN;
                medIntraCorrB(rowIdx,1:4,1:2)  = NaN;
            continue; 
        end 

       % Take average reference 
       probeARef = probeA - mean(probeA(:,chA(1):chA(2)),2,'omitnan'); 
       probeBRef = probeB - mean(probeB(:,chB(1):chB(2)),2,'omitnan'); 

       % Get within probe pairwise correlations
       corrA = max(imgaussfilt(corr(probeA),1),0);
       corrB = max(imgaussfilt(corr(probeB),1),0);

       lowIntraCorr(rowIdx,1) = mean(corrA,'all')<0.29|  mean(corrB,'all')<0.29;

        % Get mean, median and maximum pairwise correlations for different
        % frequency bands...
        clear timeSeriesCorr sizeSpecA sizeSpecB specA_R specB_R freqSpecCorr
        for iBand = 1:4
            clear xA yA xARef xBRef specA specB specEEG timeValsSpec freqValsSpec eegGood envelopeABandLimited envelopeBBandLimited infraSlowA infraSlowB
            
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||...
                    (strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))...
                    ||(strcmp(expDate,'09_19_2022') && (fileNum == 4)) || (isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate}))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||

                meanPairCorr(rowIdx,iBand,1:2)   = NaN; maxPairCorr(rowIdx,iBand,1:2)    = NaN; medPairCorr(rowIdx,iBand,1:2)  = NaN;
                meanIntraCorrA(rowIdx,iBand,1:2) = NaN; meanIntraCorrB(rowIdx,iBand,1:2) = NaN;
                medIntraCorrA(rowIdx,iBand,1:2)  = NaN; medIntraCorrB(rowIdx,iBand,1:2)  = NaN;
                continue;
            end

            if chA(1) == 0 || chB(1)== 0 % Check if any of the recordings are empty
                meanPairCorr(rowIdx,iBand,1:2)   = NaN; maxPairCorr(rowIdx,iBand,1:2)    = NaN;
                medPairCorr(rowIdx,iBand,1:2)    = NaN; meanIntraCorrA(rowIdx,iBand,1:2) = NaN;
                meanIntraCorrB(rowIdx,iBand,1:2) = NaN; medIntraCorrA(rowIdx,iBand,1:2)  = NaN;
                medIntraCorrB(rowIdx,iBand,1:2)  = NaN;
                continue;
            end

            if ~isempty(eeg{fileNum,iDate})
                eegGood                             = double(eeg{fileNum,iDate});
                eegGood(allBadTimes{fileNum,iDate}) = [];
                [eegPSD,psdFreq]                    = mtspectrumc(eegGood,params);
            end
            for iRef = 1:2
                clear xATemp yATemp 
                switch iRef
                    case 1
                        xATemp = probeA;
                        yATemp = probeB;
                    case 2
                        xATemp = probeARef;
                        yATemp = probeBRef;
                end

                % Get spectrograms for LFP (channels inside cortex) and EEG
                [specA,timeValsSpec,freqValsSpec] = mtspecgramc(xATemp,[5 2],params);
                [specB,~,~]                       = mtspecgramc(yATemp,[5 2],params);

                switch iBand % Get the correlation vs connectivity vs distance plots for different bands...
                    case 1 % Alpha band
                        xA = filtfilt(sosA,gA,double(xATemp(:,chA(1):chA(2))));
                        yA = filtfilt(sosA,gA,double(yATemp(:,chB(1):chB(2))));
                       
                        fInd  = freqValsSpec>=8 & freqValsSpec<=12; % Alpha

                        if ~isempty(eeg{fileNum,iDate})
                            fIndEEG = psdFreq>=8 & psdFreq<=12;
                        end

                    case 2 % Beta band
                        xA = filtfilt(sosB,gB,double(xATemp(:,chA(1):chA(2))));
                        yA = filtfilt(sosB,gB,double(yATemp(:,chB(1):chB(2))));

                        fInd = freqValsSpec>=13 & freqValsSpec<=30; % Beta

                        if ~isempty(eeg{fileNum,iDate})
                            fIndEEG = psdFreq>=13 & psdFreq<=30;
                        end

                    case 3% Gamma band
                        xA = filtfilt(sosG,gG,double(xATemp(:,chA(1):chA(2))));
                        yA = filtfilt(sosG,gG,double(yATemp(:,chB(1):chB(2))));

                        fInd = freqValsSpec>=30 & freqValsSpec<=90; % Gamma

                        if ~isempty(eeg{fileNum,iDate})
                            fIndEEG = psdFreq>=30 & psdFreq<= 90;
                        end

                    case 4 % Wideband
                        xA = xATemp(:,chA(1):chA(2));
                        yA = yATemp(:,chB(1):chB(2));

                        fInd = true(1,length(freqValsSpec)); % Wideband

                        if ~isempty(eeg{fileNum,iDate})
                            fIndEEG = true(1,length(psdFreq));
                        end
                end

                numChA = size(xA,2);
                numChB = size(yA,2);

                % Get the intra probe correlations for channels inside the
                % cortex...
                meanIntraCorrA(rowIdx,iBand,iRef) = mean(corr(xA,'rows','complete'),'all','omitnan');
                meanIntraCorrB(rowIdx,iBand,iRef) = mean(corr(yA,'rows','complete'),'all','omitnan');
                medIntraCorrA(rowIdx,iBand,iRef)  = median(corr(xA,'rows','complete'),'all','omitnan');
                medIntraCorrB(rowIdx,iBand,iRef)  = median(corr(yA,'rows','complete'),'all','omitnan');


                % Get pairwise correlations between the two probes...
                maxPairCorr(rowIdx,iBand,iRef)  = max(corr(xA,yA),[],'all','omitnan');
                meanPairCorr(rowIdx,iBand,iRef) = mean(corr(xA,yA),'all','omitnan');
                medPairCorr(rowIdx,iBand,iRef)  = median(corr(xA,yA),'all','omitnan');

                if isempty(timeValsSpec); continue; end

                % Get spectrogram and LFP powers for different probes for all
                % channels and get the mean/median powers across time and
                % frequencies:
                meanSpecA{fileNum,iDate}(:,iBand,iRef) = squeeze(mean(10.*log10(abs(specA(:,fInd,:))),[1 2],'omitnan'));
                meanSpecB{fileNum,iDate}(:,iBand,iRef) = squeeze(mean(10.*log10(abs(specB(:,fInd,:))),[1 2],'omitnan'));
                medSpecA{fileNum,iDate}(:,iBand,iRef)  = squeeze(median(10.*log10(abs(specA(:,fInd,:))),[1 2],'omitnan'));
                medSpecB{fileNum,iDate}(:,iBand,iRef)  = squeeze(median(10.*log10(abs(specB(:,fInd,:))),[1 2],'omitnan'));

                % Get overall powers from spectrogram for channels inside
                % cortex...
                specAMeanAll(rowIdx,iBand,iRef) = mean(meanSpecA{fileNum,iDate}(chA(1):chA(2),iBand,iRef),'omitnan');
                specBMeanAll(rowIdx,iBand,iRef) = mean(meanSpecB{fileNum,iDate}(chB(1):chB(2),iBand,iRef),'omitnan');
                specAMedAll(rowIdx,iBand,iRef)  = mean(medSpecA{fileNum,iDate}(chA(1):chA(2),iBand,iRef),'omitnan');
                specBMedAll(rowIdx,iBand,iRef)  = mean(medSpecB{fileNum,iDate}(chB(1):chB(2),iBand,iRef),'omitnan');

                % Get instantaneous power and correlate the powers
                envelopeABandLimited = envelope(abs(xA));
                envelopeBBandLimited = envelope(abs(yA));  

                % Correlate instantaneous band power
                meanCorrEnvelope(rowIdx,iBand,iRef) =  mean(corr(envelopeABandLimited,envelopeBBandLimited),'all','omitnan');
                medCorrEnvelope(rowIdx,iBand,iRef)  =  median(corr(envelopeABandLimited,envelopeBBandLimited),'all','omitnan');

                % Filter envelope from 0.01-0.1 Hz
                envelopeABandLimited = [envelopeABandLimited(1:1e3,:); envelopeABandLimited ;envelopeABandLimited(1:1e3,:) ];
                infraSlowA = filtfilt(sos,g,double(envelopeABandLimited));
                infraSlowA = single(infraSlowA(1e3+1:(end-1e3),:));

                envelopeBBandLimited = [envelopeBBandLimited(1:1e3,:); envelopeBBandLimited ;envelopeBBandLimited(1:1e3,:)];
                infraSlowB = filtfilt(sos,g,double(envelopeBBandLimited));
                infraSlowB = single(infraSlowB(1e3+1:(end-1e3),:));                

                % Correlate infraslow flucutuations in instantaneous band power
                meanCorrInfraSlow(rowIdx,iBand,iRef) =  mean(corr(infraSlowA,infraSlowB),'all','omitnan');
                medCorrInfraSlow(rowIdx,iBand,iRef)  =  median(corr(infraSlowA,infraSlowB),'all','omitnan');

                % Marginal correlations of A relative to channel B
                marginalA_B_temp{rowIdx}(:,iBand,iRef) = mean(corr(infraSlowA,infraSlowB),2,'omitnan');
                marginalB_A_temp{rowIdx}(:,iBand,iRef)  = mean(corr(infraSlowB,infraSlowA),2,'omitnan');

                % One-one channel mapping
                corrTemp = corr(infraSlowA,infraSlowB);
%                 one_oneCorr_temp{rowIdx}(iBand,iRef) = zeros(1, size(corrTemp,2));
                X = 1: min(size(corrTemp)); Y = 1:min(size(corrTemp));
                for k = 1 :min(size(corrTemp))
                    row = Y(k);
                    col = X(k);
                    one_oneCorr_temp{rowIdx}(k,iBand,iRef) = corrTemp(row, col);
                end

                % Correlate superficial and deep channels separately...
                if mod(numChA,2) == 0; sChA = floor(numChA/2);else sChA = floor((numChA+1)/2); end
                if mod(numChB,2) == 0; sChB = floor(numChB/2);else sChB = floor((numChB+1)/2); end

                superMeanCorr(rowIdx,iBand,iRef) = mean(corr(xA(:,1:sChA),yA(:,1:sChB)),'all','omitnan');
                superMedCorr(rowIdx,iBand,iRef)  = median(corr(xA(:,1:sChA),yA(:,1:sChB)),'all','omitnan');

                superMeanEnvelopeCorr(rowIdx,iBand,iRef) = mean(corr(envelopeABandLimited(:,1:sChA),envelopeBBandLimited(:,1:sChB)),'all','omitnan');
                superMedEnvelopeCorr(rowIdx,iBand,iRef)  = median(corr(envelopeABandLimited(:,1:sChA),envelopeBBandLimited(:,1:sChB)),'all','omitnan');

                superMeanInfraSlowCorr(rowIdx,iBand,iRef) = mean(corr(infraSlowA(:,1:sChA),infraSlowB(:,1:sChB)),'all','omitnan');
                superMedInfraSlowCorr(rowIdx,iBand,iRef)  = median(corr(infraSlowA(:,1:sChA),infraSlowB(:,1:sChB)),'all','omitnan');

                if ~isempty(median(corr(xA(:,sChA+1:end),yA(:,sChB+1:end)),'all','omitnan'))
                    deepMedCorr(rowIdx,iBand,iRef)  = median(corr(xA(:,sChA+1:end),yA(:,sChB+1:end)),'all','omitnan');
                    deepMeanCorr(rowIdx,iBand,iRef) = mean(corr(xA(:,sChA+1:end),yA(:,sChB+1:end)),'all','omitnan');
                    deepMeanEnvelopeCorr(rowIdx,iBand,iRef) = mean(corr(envelopeABandLimited(:,sChA+1:end),envelopeBBandLimited(:,sChB+1:end)),'all','omitnan');
                    deepMedEnvelopeCorr(rowIdx,iBand,iRef)  = median(corr(envelopeABandLimited(:,sChA+1:end),envelopeBBandLimited(:,sChB+1:end)),'all','omitnan');
                    deepMeanInfraSlowCorr(rowIdx,iBand,iRef) = mean(corr(infraSlowA(:,sChA+1:end),infraSlowB(:,sChB+1:end)),'all','omitnan');
                    deepMedInfraSlowCorr(rowIdx,iBand,iRef)  = median(corr(infraSlowA(:,sChA+1:end),infraSlowB(:,sChB+1:end)),'all','omitnan');

                else
                    deepMedCorr(rowIdx,iBand,iRef)  = superMeanCorr(rowIdx,iBand,iRef);
                    deepMeanCorr(rowIdx,iBand,iRef) = superMedCorr(rowIdx,iBand,iRef);
                    deepMeanEnvelopeCorr(rowIdx,iBand,iRef) = superMeanEnvelopeCorr(rowIdx,iBand,iRef);
                    deepMedEnvelopeCorr(rowIdx,iBand,iRef)  = superMedEnvelopeCorr(rowIdx,iBand,iRef);
                    deepMeanInfraSlowCorr(rowIdx,iBand,iRef) = superMeanInfraSlowCorr(rowIdx,iBand,iRef);
                    deepMedInfraSlowCorr(rowIdx,iBand,iRef)  = superMedInfraSlowCorr(rowIdx,iBand,iRef);
                end

                if exist('eegGood','var') && iRef == 1
                    meanPSDEEG(rowIdx,iBand) = squeeze(mean(10.*log10(abs(eegPSD(fIndEEG))),'omitnan'));
                    medPSDEEG(rowIdx,iBand)  = squeeze(median(10.*log10(abs(eegPSD(fIndEEG))),'omitnan'));
                end
            end
        end
        rowIdx = rowIdx+1;
    end
end

toc;
clear probeA probeB xA yA  probeBBandLimited probeABandLimited
% Remove bad data points from matrices...
nanIdx = find(isnan(meanPairCorr(:,1)));

% connValsAll(nanIdx)=[];
% distSitesAll(nanIdx) = [];
% heartRateValsAll(nanIdx) = [];
% anesthesiaValsAll(nanIdx) = [];
% meanIntraCorrA(nanIdx,:,:)=[];
% meanIntraCorrB(nanIdx,:,:)=[];
% medIntraCorrA(nanIdx,:,:)=[];
% medIntraCorrB(nanIdx,:,:)=[];
% 
% meanPairCorr(nanIdx,:,:) = [];
% medPairCorr(nanIdx,:,:)  = [];
% maxPairCorr(nanIdx,:,:)  = [];
% 
% specAMeanAll(nanIdx,:,:) = []; specBMeanAll(nanIdx,:,:) = [];
% specAMedAll(nanIdx,:,:)  = []; specBMedAll(nanIdx,:,:) = [];
% meanPSDEEG(nanIdx,:)  = []; medPSDEEG(nanIdx,:) = [];
% 
% meanCorrEnvelope(nanIdx,:,:) = []; medCorrEnvelope(nanIdx,:,:) = [];
% meanCorrInfraSlow(nanIdx,:,:) = []; medCorrInfraSlow(nanIdx,:,:) = [];
% % meanCorrInfraSlowLFP(nanIdx,:) = []; medCorrInfraSlowLFP(nanIdx,:) = []; 
% 
% superMeanCorr(nanIdx,:,:) = []; superMedCorr(nanIdx,:,:)= []; 
% superMeanEnvelopeCorr(nanIdx,:,:) = []; superMedEnvelopeCorr(nanIdx,:,:)= [];
% superMeanInfraSlowCorr(nanIdx,:,:) = []; superMedInfraSlowCorr(nanIdx,:,:) = []; 
% 
%  deepMedCorr(nanIdx,:,:)  = [];  deepMeanCorr(nanIdx,:) = [];
%  deepMeanEnvelopeCorr(nanIdx,:,:) = [];  deepMedEnvelopeCorr(nanIdx,:,:)  = [];
%  deepMeanInfraSlowCorr(nanIdx,:,:) =[]; deepMedInfraSlowCorr(nanIdx,:,:)  = [];

% movMeanPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:) = [];
save(['D:\Data\' monkeyName '_SqM\Left Hemisphere\' monkeyName '_varsNew.mat'],'allBadTimes','badElecA','badElecB','estChInCortexA','estChInCortexB',...
    'connValsAll','distSitesAll','meanIntraCorrA','meanIntraCorrB', 'medIntraCorrA', 'medIntraCorrB','meanPairCorr','medPairCorr','specAMeanAll','specBMeanAll',...
    'meanPSDEEG','medPSDEEG','meanCorrEnvelope','medCorrEnvelope','meanCorrInfraSlow','medCorrInfraSlow','superMeanCorr','superMedCorr', 'superMeanEnvelopeCorr',...
    'superMedEnvelopeCorr','superMeanInfraSlowCorr', 'superMedInfraSlowCorr','deepMedCorr','deepMeanCorr','deepMeanEnvelopeCorr','deepMedEnvelopeCorr',...
    'deepMeanInfraSlowCorr','deepMedInfraSlowCorr','marginalA_B_temp','marginalB_A_temp','one_oneCorr_temp','pairClass','patchInfo','nanIdx');











