%% ephysDualProbeRS_v3
% This program performs analysis on simultaneously recorded LFP data from
% two laminar electrode arrays. 
% March 15,2023 - Keerthana Manikandan
% Refer ephysDualProbeRS_v1 and ephysDualProbeRS_V2 for previous iterations
% of analyses performed. This program performs the following - 
% 1. Stores and retrieves LFP for a monkey
% 2. Removes bad time segments from the data 
% 3. Computes intra probe and pairwise correlations between the two probes
% 4. Computes the distance and connectivity between two probes. 
% 5. Determines the transition channel 
% 6. Computes the mean and median pairwise correlations for different
% frequency bands for channels inside cortex. 
% 7. Plots the relationship between pairwise correlations and functional
% connectivity for different bands 
% 8. Plots the relationship between pairwise correlations and distance for
% different bands. 
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
monkeyName     = 'Whiskey';
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
    chInCortexProbeAUpdated,chInCortexProbeBUpdated,probeLabelA,probeLabelB,anesthesiaLevels,heartRate] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']); 

%% 2. Store/Retrieve the LFP (1-250 Hz) and also removing 1-5 Hz
% probe1 = {}; probe2 ={}; eeg = {}; scalpEEG = {};
disp(['Storing/Retrieving LFP data for ' monkeyName]);
for iDate = 1:size(allDates,1)
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

            if exist('eegCh','var')
                eeg{fileNum,iDate} = eegCh;
            else
                eeg{fileNum,iDate} = [];
            end

            if exist('scalpEEGCh','var')
                scalpEEG{fileNum,iDate} = scalpEEGCh;
            else
                scalpEEG{fileNum,iDate} = [];
            end
        end
    end
end
clc; disp('Data Stored/Retrieved');

%% 3.Get/retrieve the distance and the connectivity values for the sites
clear distSites connSites greenMapRef 
[distSites,connVals,refSites,movSites,greenMapRef] = getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);

% Get the connectivity and distance in a vector
distSitesAll = []; connValsAll = []; heartRateValsAll = []; anesthesiaValsAll = [];

for iDate = 1:size(allDates,1) 
    datFileNum        = datFileNumAll{iDate,1};
    distSitesAll      = [distSitesAll;distSites{iDate,1}(datFileNum)];
    connValsAll       = [connValsAll; squeeze(mean(connVals{iDate,1}(:,:,datFileNum),[1,2],'omitnan'))];
    anesthesiaValsAll = [anesthesiaValsAll; anesthesiaLevels{iDate,1}(datFileNum)]; 
    heartRateValsAll  = [heartRateValsAll; heartRate{iDate,1}(datFileNum)];
end
disp(['Obtained/retrieved distance between probes and connectivity values for ' monkeyName]);

%% 4. Remove bad time segments from data and store intra probe, pairwise correlograms
saveFigureFlag = 0;
%  [removeTimes] = getCorrelograms(monkeyName, allDates, datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,saveFigureFlag);
[allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);
% [probe1,probe2,badElecA,badElecB] = getCorrelograms_vTemp(monkeyName, allDates, datFileNumAll,probe1,probe2,probeLabelA,probeLabelB,saveFigureFlag);

%% 5. Get the transition channels for the recordings... 
% Step 1: Check powers to see if there are bad channels.
% Step 2: Eliminate bad channels and obtain the transition channels by
% computing slope of the marginals (obtained by averaging the intra-probe
% correlograms) - Mean is used here for sensitivity to outliers which is
% needed for picking the transition channel

for iDate = 1:size(allDates,1) % all experiment dates
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
        
        % Remove extra channels if present
        if size(probe1{fileNum,iDate},2) == 33; probe1{fileNum,iDate}(:,1) = []; end
        if size(probe2{fileNum,iDate},2) == 33; probe2{fileNum,iDate}(:,1) = []; end


      % Obtain data for Probe A and Probe B
        probeA = probe1{fileNum,iDate};%filtfilt(bG,aG,probe1{fileNum,iDate});
        probeB = probe2{fileNum,iDate}; %filtfilt(bG,aG,probe2{fileNum,iDate}); 

        % Remove bad channels from Probe B
         if chInCortexProbeB{iDate}(iFile)~= 1
            probeB(:,badChProbeB)= [];
        end

        % Remove bad channels
        if ~isempty(badElecA{fileNum,iDate})
            probeA(:,badElecA{fileNum,iDate}) = [];
        end

        if ~isempty(badElecB{fileNum,iDate})
            probeB(:,ismember(probeBList,badElecB{fileNum,iDate})) = [];
        end  

        % Remove the bad time segments
        if ~isempty(allBadTimes{fileNum,iDate})
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];
        end
      
        % Step 2: Get the intra probe marginals from wideband range
        marginalA = mean(corr(probeA,'Rows','complete'),2,'omitnan');

        fxA = abs(movmean(gradient(marginalA),2,'omitnan')); % Get the slope of the marginals
        if strcmp(expDate,'11_01_2021') || strcmp(expDate,'01_01_2021')
            estChInCortexA{iDate}(fileNum,1) = 1;
        else
            estChInCortexA{iDate}(fileNum,1) = find(fxA == max(fxA)); % Find the channel with maximum slope
        end

        if estChInCortexA{iDate}(fileNum,1)==20
            estChInCortexA{iDate}(fileNum,1) = (estChInCortexA{iDate}(fileNum,1) - 19);
        elseif estChInCortexA{iDate}(fileNum,1)>20
            estChInCortexA{iDate}(fileNum,1) = (estChInCortexA{iDate}(fileNum,1) - 20);
        end

        estChInCortexA{iDate}(fileNum,2) = size(probeA,2);

        % Repeat the above process for Probe B
        marginalB = mean(corr(probeB,'Rows','complete'),2,'omitnan');
        fxB = abs(movmean(gradient(marginalB),2,'omitnan'));

        estChInCortexB{iDate}(fileNum,1) = find(fxB == max(fxB));

        if estChInCortexB{iDate}(fileNum,1) == 20
            estChInCortexB{iDate}(fileNum,1) = (estChInCortexB{iDate}(fileNum,1) - 19);
        elseif estChInCortexB{iDate}(fileNum,1)>20
            estChInCortexB{iDate}(fileNum,1) = (estChInCortexB{iDate}(fileNum,1) - 20);
        end

        estChInCortexB{iDate}(fileNum,2) = size(probeB,2);

%         
%         % Method 2: Perform the same steps except taking the gamma band
%         % into consideration 
%         clear marginalA fxA marginalB fxB 
%         marginalA = mean(corr(filtfilt(bG,aG,probeA),'Rows','complete'),2,'omitnan');
% 
%         fxA = abs(movmean(gradient(marginalA),2,'omitnan')); % Get the slope of the marginals
%         estChInCortexAGamma{iDate}(fileNum,1) = find(fxA == max(fxA))+1; % Find the channel with maximum slope
% 
%         if estChInCortexAGamma{iDate}(fileNum,1)>=20
%             estChInCortexAGamma{iDate}(fileNum,1) = (estChInCortexAGamma{iDate}(fileNum,1) - 19);
%         end
% 
%         if estChInCortexAGamma{iDate}(fileNum,1)+19>size(probeA,2) || strcmp(expDate,'10_17_2022') % Probe A on this date has spacing 50um
%             estChInCortexAGamma{iDate}(fileNum,2) = size(probeA,2);
%         else
%             estChInCortexAGamma{iDate}(fileNum,2) = estChInCortexAGamma{iDate}(fileNum,1)+19;
%         end
% 
%         % Repeat the above process for Probe B
%         marginalB = mean(corr(filtfilt(bG,aG,probeB),'Rows','complete'),2,'omitnan');
%         fxB = abs(movmean(gradient(marginalB),2,'omitnan'));
% 
%         estChInCortexBGamma{iDate}(fileNum,1) = find(fxB == max(fxB))+1;
% 
%         if estChInCortexBGamma{iDate}(fileNum,1)>=20
%             estChInCortexBGamma{iDate}(fileNum,1) = (estChInCortexBGamma{iDate}(fileNum,1) - 19);
%         end
% 
%         if estChInCortexBGamma{iDate}(fileNum,1)+19>size(probeB,2)
%             estChInCortexBGamma{iDate}(fileNum,2) = size(probeB,2);
%         else
%             estChInCortexBGamma{iDate}(fileNum,2) = estChInCortexBGamma{iDate}(fileNum,1)+19;
%         end

        % Method 3: Get the powers in 40-120 Hz range and look at median power for
        % each channel 
%         clear freqValsSpec spec medHighFreqSpec fxG 
%         [spec,~,freqValsSpec] = mtspecgramc(probeA,[5 2],params);
%         fInd = freqValsSpec>=40;
%         medHighFreqSpec = squeeze(median(10.*log10(abs(spec(:,fInd,:))),[1 2],'omitnan'));
%         fxG = abs(movmean(gradient(medHighFreqSpec),2,'omitnan'));
%         estChInCortexASpec{iDate}(fileNum,1) = find(abs(fxG) == max(abs(fxG)))+1;
% 
%         clear freqValsSpec spec medHighFreqSpec fxG
%         [spec,~,freqValsSpec] = mtspecgramc(probeB,[5 2],params);
%         fInd = find(freqValsSpec>=40);
%         medHighFreqSpec = squeeze(median(10.*log10(abs(spec(:,fInd,:))),[1 2],'omitnan'));
%         fxG = abs(movmean(gradient(medHighFreqSpec),2,'omitnan'));
%         estChInCortexBSpec{iDate}(fileNum,1) = find(abs(fxG) == max(abs(fxG)))+1;

    end
end 

%% 6. Get the powers - instantaneous power, bandlimited power, spectrogram and pairwise correlations for each animal 
% clear meanPairCorr medPairCorr maxPairCorr meanIntraCorrA meanIntraCorrB medIntraCorrA medIntraCorrB
rowIdx = 1; eegFlag = []; 
bandLabels = {'Alpha band'; 'Beta band'; 'Gamma band';'Wideband'}; 
tic;
for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum)
%         if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
%         if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end
        fileNum = datFileNum(iFile);

        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
     
        if ~isempty(eeg{fileNum,iDate}); eegFlag(rowIdx) = 1; else eegFlag(rowIdx) = 0; end 

        if isempty(probeA) || isempty(probeB)
            meanPairCorr(rowIdx,:) = NaN;
            maxPairCorr(rowIdx,:)  = NaN;
            medPairCorr(rowIdx,:)  = NaN;

            meanIntraCorrA(rowIdx,:) = NaN;
            meanIntraCorrB(rowIdx,:) = NaN;
            medIntraCorrA(rowIdx,:)  = NaN;
            medIntraCorrB(rowIdx,:)  = NaN;

            rowIdx = rowIdx+1;
            continue;
        end

        if size(probeB,2) == 33; probeB(:,1) = []; end
        if size(probeA,2) == 33; probeA(:,1) = []; end
        

           % Remove bad channels from Probe B
         if chInCortexProbeB{iDate}(iFile)~= 1
            probeB(:,badChProbeB)= [];
        end

        % Remove channels that have high impedance - from datasheets
        if ~isempty(badElecA{fileNum,iDate})
            probeA(:,badElecA{fileNum,iDate}) = [];
        end
        if ~isempty(badElecB{fileNum,iDate})
            probeB(:,ismember(probeBList,badElecB{fileNum,iDate})) = [];
        end  

        % Remove the bad time segments
        if ~isempty(allBadTimes{fileNum,iDate})
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];
        end


        % Get power - instantaneous, power spectrum and spectrogram
%         [powA{fileNum,iDate},powB{fileNum,iDate},powEEG{fileNum,iDate},specA{fileNum,iDate},specB{fileNum,iDate}...
%             ,specEEG{fileNum,iDate},instPowA{fileNum,iDate},instPowB{fileNum,iDate},freqVals,freqValsSpec,timeValsSpec] =...
%             getAllPowerData(probeA,probeB,eeg{fileNum,iDate},fs,monkeyName,expDate,fileNum);

        % Get pairwise correlations
        corrA = max(imgaussfilt(corr(probeA),1),0);
        corrB = max(imgaussfilt(corr(probeB),1),0);
        lowIntraCorr(rowIdx,1) = mean(corrA,'all')<0.29|  mean(corrB,'all')<0.29;
        
%         disp('Obtaining pairwise correlations...');
        % Get mean, median and maximum pairwise correlations for different
        % frequency bands... 
        for iBand = 1:4
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||(strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))||(strcmp(expDate,'09_19_2022') && (fileNum == 4))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||
                meanPairCorr(rowIdx,iBand) = NaN;
                maxPairCorr(rowIdx,iBand)  = NaN;
                medPairCorr(rowIdx,iBand)  = NaN;

                meanIntraCorrA(rowIdx,iBand) = NaN;
                meanIntraCorrB(rowIdx,iBand) = NaN;
                medIntraCorrA(rowIdx,iBand)  = NaN;
                medIntraCorrB(rowIdx,iBand)  = NaN;
                continue;
            end

            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate})
                meanPairCorr(rowIdx,iBand) = NaN; maxPairCorr(rowIdx,iBand) = NaN;
                medPairCorr(rowIdx,iBand)  = NaN; 
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
%             if chA(1)+20>size(probeA,2)
                chA(2) = size(probeA,2);
%             else
%                 chA(2) = chA(1)+20; 
%             end 

            chB(1) = estChInCortexB{iDate}(fileNum,1);

            if chB(1) == 0
                meanPairCorr(rowIdx,iBand) = NaN;maxPairCorr(rowIdx,iBand) = NaN;
                medPairCorr(rowIdx,iBand) = NaN; 
                continue;
            end

            if ~strcmp(expDate,'08_08_2022')
%                 if chB(1)+20>size(probeB,2)
                    chB(2) = size(probeB,2);
%                 else
%                     chB(2) = chB(1)+20;
%                 end
            else
                chB(2) = chB(1);
            end
            switch iBand % Get the correlation vs connectivity vs distance plots for different bands...
                case 4 % Wideband
                    xA = probeA(:,chA(1):chA(2));
                    yA = probeB(:,chB(1):chB(2));

                case 1 % Alpha band
                    xA = filtfilt(bA,aA,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bA,aA,probeB(:,chB(1):chB(2)));

                case 2 % Beta band
                    xA = filtfilt(bB,aB,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bB,aB,probeB(:,chB(1):chB(2)));

                case 3% Gamma band
                    xA = filtfilt(bG,aG,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bG,aG,probeB(:,chB(1):chB(2)));
            end

            % Get the intra probe correlations for channels inside the
            % cortex...
            meanIntraCorrA(rowIdx,iBand) = mean(corr(xA,'rows','complete'),'all','omitnan');
            meanIntraCorrB(rowIdx,iBand) = mean(corr(yA,'rows','complete'),'all','omitnan');
            medIntraCorrA(rowIdx,iBand)  = median(corr(xA,'rows','complete'),'all','omitnan');
            medIntraCorrB(rowIdx,iBand)  = median(corr(yA,'rows','complete'),'all','omitnan');

            % Get pairwise correlations between the two probes...
            maxPairCorr(rowIdx,iBand)  =  max(corr(xA,yA),[],'all','omitnan');
            meanPairCorr(rowIdx,iBand) = mean(corr(xA,yA),'all','omitnan');
            medPairCorr(rowIdx,iBand)  = median(corr(xA,yA),'all','omitnan');
        end 
         rowIdx = rowIdx+1;
    end
end
toc;

%%
connValsOld = connValsAll; 
distOld     = distSitesAll; 
%% Run this once only...
oneIdx = [];
nanIdx = find(isnan(meanPairCorr(:,1)));
lowIntraCorr = [];%find(lowIntraCorr==1);
% oneIdx = find(meanIntraCorrB(:,1) == 1);
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
% movMeanPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:) = [];

%%
figure;
cVal = {'b';'r';'m';'k'};
for iBand = 1:4
    subplot(2,2,iBand);
    scatter(connValsAll,medPairCorr(:,iBand),35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand});hold on; 
%     if ~isempty(~eegFlag)
%         scatter(connValsAll(~eegFlag),medPairCorr(~eegFlag,iBand),35,'filled','MarkerFaceColor','r','MarkerEdgeColor',cVal{iBand});
%     end 

    xlim([-0.3 1]); ylim([-1 1]);
    xlabel('Functional Connectivity'); ylabel('Pairwise Correlations');
    hold on; box on;
    title(bandLabels{iBand})

    xFit = linspace(min(connValsAll),max(connValsAll),1000);
    pairWiseCorr = medPairCorr(:,iBand);
    pairWiseCorr(isnan(pairWiseCorr)) = 0;
    [f,gof]  = fit(connValsAll,pairWiseCorr,'poly1');
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    r2Final (iBand) = r2; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);
end 
%%
figure;
cVal = {'b';'r';'m';'k'};
for iBand = 1:4
    subplot(2,2,iBand);
    scatter(distSitesAll(1:45),medPairCorr(1:45,iBand),35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand}); hold on
%      if ~isempty(~eegFlag)
%         scatter(distSitesAll(~eegFlag),medPairCorr(~eegFlag,iBand),35,'filled','MarkerFaceColor','r','MarkerEdgeColor',cVal{iBand});
%     end
    xlim([0 20]); ylim([-1 1]);
    xlabel('Distance (mm)'); ylabel('Pairwise Correlations');
    hold on; box on;
    title(bandLabels{iBand})

    xFit = linspace(min(distSitesAll),max(distSitesAll),1000);
    pairWiseCorr = medPairCorr(:,iBand);
    pairWiseCorr(isnan(pairWiseCorr)) = 0;
    [f,gof]  = fit(distSitesAll,pairWiseCorr,'poly2');
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit.^2) +c(2)*(xFit)+c(3);
    r2Final (iBand) = r2; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);
end 

%% Mean and median powers for different frequencies from spectrogram 
tic;
clear meanSpecEEG medSpecEEG  specAMeanAll specBMeanAll specAMedAll specBMedAll
rowIdx          = 1; 
fs              = 1e3; % Sampling frequency
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;
gammaBand       = [30 90]; 
alphaBand       = [8 12]; 
fileCount       = 1; 
for iDate = 1:size(allDates,1)
    clc; clear expdate datFileNum datFileName
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        clear probeA probeB 
        fileNum = datFileNum(iFile);
        clc; disp(['Obtaining and processing spectrograms for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);
        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); rowIdx = rowIdx+1; continue; end % Check if probes do not have any LFP
        
        % Remove extra channels if present
        if size(probe1{fileNum,iDate},2) == 33; probe1{fileNum,iDate}(:,1) = []; end
        if size(probe2{fileNum,iDate},2) == 33; probe2{fileNum,iDate}(:,1) = []; end


      % Obtain data for Probe A and Probe B
        probeA = probe1{fileNum,iDate};%filtfilt(bG,aG,probe1{fileNum,iDate});
        probeB = probe2{fileNum,iDate}; %filtfilt(bG,aG,probe2{fileNum,iDate}); 

        % Remove bad channels from Probe B
         if chInCortexProbeB{iDate}(iFile)~= 1
            probeB(:,badChProbeB)= [];
        end

        % Remove bad channels
        if ~isempty(badElecA{fileNum,iDate})
            probeA(:,badElecA{fileNum,iDate}) = [];
        end

        if ~isempty(badElecB{fileNum,iDate})
            probeB(:,ismember(probeBList,badElecB{fileNum,iDate})) = [];
        end  

        % Remove the bad time segments
        if ~isempty(allBadTimes{fileNum,iDate})
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];
        end
        
        % Get spectrograms
        clear specA specB specEEG timeValsSpec freqValsSpec eegGood 
        [specA,timeValsSpec,freqValsSpec] = mtspecgramc(probeA,[5 2],params);
        [specB,~,~] = mtspecgramc(probeB,[5 2],params);
        if ~isempty(eeg{fileNum,iDate})
            eegGood = eeg{fileNum,iDate};
            eegGood(allBadTimes{fileNum,iDate}) = []; 
            [specEEG,~,~] = mtspecgramc(eegGood,[5 2],params);           
        end 
        if isempty(timeValsSpec); rowIdx = rowIdx+1; continue; end 

        for iBand = 1:4
            clear fInd; 
            switch iBand
                case 1
                    fInd = freqValsSpec>=8 & freqValsSpec<=12; % Alpha
                case 2
                    fInd = freqValsSpec>=13 & freqValsSpec<=30; % Beta
                case 3
                    fInd = freqValsSpec>=30 & freqValsSpec<=90; % Gamma
                case 4
                    fInd = true(1,length(freqValsSpec)); % Wideband
            end
                % Get the mean/median powers across time and frequencies -
                % channelwise 
                meanSpecA{fileNum,iDate}(:,iBand) = squeeze(mean(10.*log10(abs(specA(:,fInd,:))),[1 2],'omitnan'));
                meanSpecB{fileNum,iDate}(:,iBand) = squeeze(mean(10.*log10(abs(specB(:,fInd,:))),[1 2],'omitnan'));
                medSpecA{fileNum,iDate}(:,iBand)  = squeeze(median(10.*log10(abs(specA(:,fInd,:))),[1 2],'omitnan'));
                medSpecB{fileNum,iDate}(:,iBand)  = squeeze(median(10.*log10(abs(specB(:,fInd,:))),[1 2],'omitnan'));
                
                % Get overall powers from spectrogram 
                specAMeanAll(rowIdx,iBand) = mean(meanSpecA{fileNum,iDate}(:,iBand),'omitnan');
                specBMeanAll(rowIdx,iBand) = mean(meanSpecB{fileNum,iDate}(:,iBand),'omitnan');
                specAMedAll(rowIdx,iBand)  = mean(medSpecA{fileNum,iDate}(:,iBand),'omitnan');
                specBMedAll(rowIdx,iBand)  = mean(medSpecB{fileNum,iDate}(:,iBand),'omitnan');

                if exist('eegGood','var')
                    meanSpecEEG(rowIdx,iBand) = squeeze(mean(10.*log10(abs(specEEG(:,fInd,:))),[1 2],'omitnan'));
                    medSpecEEG(rowIdx,iBand)  = squeeze(median(10.*log10(abs(specEEG(:,fInd,:))),[1 2],'omitnan'));
                end
        end 
        rowIdx = rowIdx+1; 
    end
end
clc; disp(['Obtaining and processing spectrograms for ' monkeyName ' - Completed' ]);
toc;
%%
specAMeanAll(nanIdx,:) = []; specBMeanAll(nanIdx,:) = []; 
specAMedAll(nanIdx,:)  = []; specBMedAll(nanIdx,:) = []; 
meanSpecEEG(nanIdx,:)  = []; medSpecEEG(nanIdx,:) = []; 
%% Plot Average LFP Power vs Anesthesia, LFP Power vs Heart Rate, LFP Power vs EEG Power for different frequencies
for iPlot = 3:4
    switch iPlot
        case 1
            var = anesthesiaValsAll; 
            varName = 'Anesthesia Levels (%)';
        case 2
            var = heartRateValsAll; 
            varName = 'Heart rate (bpm)';

        case 3
            var = meanSpecEEG; 
            varName = 'Mean EEG Power (dB)';
        case 4
            var = medSpecEEG; 
            varName = 'Median EEG Power (dB)'; 
    end 

    if iPlot >=3
        if strcmp(monkeyName,'Whiskey')
             var(1:4 ,:) = [];
%             var([1:4 20:22],:) = [];
%         elseif strcmp(monkeyName,'CharlieSheen')
%             var(18:33,:) = []; 
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
                lfpA(1:4 ,:) = []; lfpB(1:4,:) = [];
                %lfpA([1:4 20:22],:) = []; lfpB([1:4 20:22],:) = [];
%             elseif strcmp(monkeyName,'CharlieSheen')
%                 lfpA(18:33 ,:) = []; lfpB(18:33,:) = [];
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

            if strcmp(monkeyName,'Whiskey')
                scatter(var(6,iBand),lfpA(6,iBand),80,'yellow','filled'); % Burst suprression
                scatter(var(7,iBand),lfpA(7,iBand),80,'green','filled');  % Deep 

                scatter(var(6,iBand),lfpB(6,iBand),80,[0.929 0.694 0.125],'filled');
                scatter(var(7,iBand),lfpB(7,iBand),80,[ 0.466 0.674 0.188],'filled');
            end 
            title(bandLabels{iBand});
            xlabel(varName); ylabel(lfpTitle); box on; grid on;
            if iBand == 1; legend('Probe A','Probe B','Location','northeast'); end
                varTemp = var; 
                aTemp = lfpA;
                bTemp = lfpB;

            if iPlot==2%>=2
              
                if iPlot == 2   
                    nanHR = find(isnan(varTemp));
                    varTemp(nanHR) = []; aTemp(nanHR,:) = []; bTemp(nanHR,:) = [];
                    xFit = linspace(min(varTemp),max(varTemp),1000);
                    pairWiseCorr = aTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fA,~]  = fit(varTemp,pairWiseCorr,'poly1');
                    pairWiseCorr =bTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fB,~]  = fit(varTemp,pairWiseCorr,'poly1');
                else
                    xFit = linspace(min(varTemp(:,iBand)),max(varTemp(:,iBand)),1000);
                    pairWiseCorr = aTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fA,~]  = fit(varTemp(:,iBand),pairWiseCorr,'poly1');
                    pairWiseCorr = bTemp(:,iBand);pairWiseCorr(isnan(pairWiseCorr)) = 0;
                    [fB,~]  = fit(varTemp(:,iBand),pairWiseCorr,'poly1');
                end
                          
                c        = coeffvalues(fA);
                fitLine   = c(1)*(xFit) +c(2);
                plot(xFit,fitLine,'Color','red','LineStyle','--','LineWidth',2);  

                c        = coeffvalues(fB);
                fitLine   = c(1)*(xFit) +c(2);
                plot(xFit,fitLine,'Color','blue','LineStyle','--','LineWidth',2);

            end 
        end
        sgtitle([lfpTitle ' vs ' varName]);
    end
end 

%% Plot boxplots for each anesthesia level for both probes A and probe B
powA = []; powB = []; groupVal = [];
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

    elseif strcmp(monkeyName,'Whiskey') && iL>4
        if anesthesiaValsAll(iL)<=1.2
            groupVal = [groupVal; 1];
            powA = [powA; specAMedAll(iL,4)];
            powB = [powB; specBMedAll(iL,4)];

        elseif any(anesthesiaValsAll(iL)>1.2 & anesthesiaValsAll(iL)<2)
            powA = [powA; specAMedAll(iL,4)];
            powB = [powB; specBMedAll(iL,4)];
            groupVal = [groupVal; 2];
          
        end
    end
end


% for iL = 1:length(uniqueAnesthesia)
%     powA{iL} = [specAMedAll((index==iL),4)];
%     powB{iL} = [specBMedAll((index==iL),4)];
% end


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

        % Remove bad channels from Probe B
        if chInCortexProbeB{iDate}(iFile)~= 1
            probeB(:,badChProbeB)= [];
        end

        % Remove bad channels
        if ~isempty(badElecA{fileNum,iDate})
            probeA(:,badElecA{fileNum,iDate}) = [];
        end

        if ~isempty(badElecB{fileNum,iDate})
            probeB(:,ismember(probeBList,badElecB{fileNum,iDate})) = [];
        end

        % Remove the bad time segments
        if ~isempty(allBadTimes{fileNum,iDate})
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];
        end

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
            infraA = filtfilt(sos,g,envelopeA);
            infraB = filtfilt(sos,g,envelopeB);

            % Get pairwise correlations
            pairCorrInfra(rowIdx,iBand) = median(corr(infraA,infraB),'all','omitnan');
        end
        rowIdx = rowIdx+1; 
    end
end 
toc;
%% Time lagged cross correlations 
rowIdx = 1; 
bandLabels = {'Alpha band'; 'Beta band'; 'Gamma band';'Wideband'}; 
tic;
for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum)

       
%         if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
%         if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end
        fileNum = datFileNum(iFile);

        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        if isempty(probeA) || isempty(probeB); continue; end 

        if isempty(probeA) || isempty(probeB)
            corrValsLagged(rowIdx,iBand) = NaN;
            lags(rowIdx,iBand) = NaN;
            rowIdx = rowIdx+1;
            continue; 
        end

        if size(probeB,2) == 33; probeB(:,1) = []; end
        if size(probeA,2) == 33; probeA(:,1) = []; end

        % Get pairwise correlations
        corrA = max(imgaussfilt(corr(probeA),1),0);
        corrB = max(imgaussfilt(corr(probeB),1),0);
        lowIntraCorr(rowIdx,1) = mean(corrA,'all')<0.29|  mean(corrB,'all')<0.29;
        
%         disp('Obtaining pairwise correlations...');
        % Get mean, median and maximum pairwise correlations for different
        % frequency bands... 
        for iBand = 1:4
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||(strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))||(strcmp(expDate,'09_19_2022') && (fileNum == 4))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||
                corrValsLagged(rowIdx,iBand) = NaN;
                lags(rowIdx,iBand) = NaN;
                continue;
            end

            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate})
                corrValsLagged(rowIdx,iBand) = NaN;
                lags(rowIdx,iBand) = NaN;
                continue;
            end

            clear xA yA chA chB
            chA = estChInCortexA{iDate}(fileNum,:);

            if chA(1) == 0
                corrValsLagged(rowIdx,iBand) = NaN;
                lags(rowIdx,iBand) = NaN;
                continue;
            end

            chB = estChInCortexB{iDate}(fileNum,:);
            if chB(1) == 0
                 corrValsLagged(rowIdx,iBand) = NaN;
                lags(rowIdx,iBand) = NaN;
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
            clear corrCh lagsCh

            for iC = 1:min([sizeA sizeB])
                clear corrValsT lagsT
                [corrValsT,lagsT] = xcorr(xA(:,iC),yA(:,iC),50,'normalized');
                corrCh(iC) = max(corrValsT);
                lagsCh(iC) = lagsT(corrValsT == max(corrValsT));  % Pick the lag where the maximum correlation occurred.
            end 

            corrValsLagged(rowIdx,iBand) = median(corrCh,'omitnan'); % Taking the median correlations
            lags(rowIdx,iBand)           = median(lagsCh,'omitnan'); % Taking the median lag 
        end 

         rowIdx = rowIdx+1;
    end
   
end
%% Do this once only

nanVals = find(isnan(corrValsLagged(:,1)));
connValsAll(nanVals)=[];
corrValsLagged(nanVals,:) =[];
lags(nanVals,:) = [];

%%
figure; 
for iBand = 1:4
  
    subplot(2,2,iBand);
    scatter(connValsAll,squeeze(corrValsLagged(:,iBand)),20,'filled');
    xlabel('Functional connectivity'); xlim([-0.3 1]);
    ylabel('Pairwise correlations'); ylim([-1 1]);
    hold on; box on;title(bandLabels{iBand}); grid on;

    xFit = linspace(min(connValsAll),max(connValsAll),1000);
    [f,gof]  = fit(connValsAll,squeeze(corrValsLagged(:,iBand)),'poly1');

    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);
    text(0.6,-0.8,['R^2: ' num2str(r2)]);
end 

f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\LaggedCorr.png'],'Resolution',300);
close gcf;



%% 
figure; 
for iBand = 1:4
    subplot(2,2,iBand);
    scatter(lags(:,iBand),squeeze(corrValsLagged(:,iBand)),10,'filled');
    xlabel('Lags'); xlim([-40 40]);
    ylabel('Pairwise correlations'); ylim([-1 1]);
    hold on; box on;title(bandLabels{iBand}); grid on;
end 
f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\Lags.png'],'Resolution',300);
close gcf;

%%
figure; scatter(distSitesAll,lags(:,4),30,'filled');
xlabel('Distance (mm)'); ylabel('Lags'); ylim([-40 40]);
box on; grid on; xlim([0 20]); 
f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\DistanceVsLag.png'],'Resolution',300);
close gcf;

%% Remove lags outside +10 and -10 
noBigLag = find(~(lags(:,4)>10 | lags(:,4)<-10));
figure; 
for iBand = 1:4
  
    subplot(2,2,iBand);
    scatter(connValsAll(noBigLag),squeeze(corrValsLagged(noBigLag,iBand)),20,'filled');
    xlabel('Functional connectivity'); xlim([-0.3 1]);
    ylabel('Pairwise correlations'); ylim([-1 1]);
    hold on; box on;title(bandLabels{iBand}); grid on;

    xFit = linspace(min(connValsAll(noBigLag)),max(connValsAll(noBigLag)),1000);
    [f,gof]  = fit(connValsAll(noBigLag),squeeze(corrValsLagged(noBigLag,iBand)),'poly1');

    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);
    text(0.6,-0.8,['R^2: ' num2str(r2)]);
end 

f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\LaggedCorrConstrained.png'],'Resolution',300);
close gcf;

%% Benchmarking the fit by shuffling Probe A/ Probe B time series
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

        if size(probeB,2) == 33; probeB(:,1) = []; end
        if size(probeA,2) == 33; probeA(:,1) = []; end

        
%         disp('Obtaining pairwise correlations...');
        % Get mean, median and maximum pairwise correlations for different
        % frequency bands... 
        for iBand = 1:4
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||(strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))||(strcmp(expDate,'09_19_2022') && (fileNum == 4))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||
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

%%
plotVals = pairCorrTimeSeries; 
colorVals =[0.6350 0.0780 0.1840];% [0 0.4470 0.7410];% This is for Whiskey
figure; 
for iBand = 1:4
    subplot(2,2,iBand);
%     for iPr = 1:2
%         clear plotVals colorVals
%         switch iPr
%             case 1
%                 plotVals = corrValsShuffledA; 
%                 colorVals = [0.4940 0.1840 0.5560];
%             case 2
%                 plotVals = corrValsShuffledB;
%                 colorVals = [0.85 0.325 0.098];
%         end

        scatter(connValsAll,squeeze(plotVals(:,iBand)),20,colorVals,'filled');
        xlabel('Functional connectivity'); xlim([-0.3 1]);
        ylabel('Pairwise correlations'); ylim([-1 1]);
        hold on; box on;title(bandLabels{iBand}); grid on;

        xFit = linspace(min(connValsAll),max(connValsAll),1000);
        [f,gof]  = fit(connValsAll,squeeze(plotVals(:,iBand)),'poly1');

        c        = coeffvalues(f);
        r2       = gof.rsquare;
        fitLine   = c(1)*(xFit) +c(2);
        plot(xFit,fitLine,'Color',colorVals,'LineStyle','--','LineWidth',2);
        text(0.6,-0.6-0.1*(iPr-1),['R^2: ' num2str(r2)],'Color',colorVals);
%     end 
%     legend('Probe A shuffled','' ,'Probe B shuffled','');
end 

% f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results\LaggedCorr.png'],'Resolution',300);
% close gcf;

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