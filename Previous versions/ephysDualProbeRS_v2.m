%% ephysDualProbeRS_v2
% This program performs analysis on simultaneously recorded LFP data from
% two laminar electrode arrays. 
% October 19,2022 - Keerthana Manikandan
% Refer ephysDualProbeRS_v1 for previous iterations of analyses performed. 
% This program performs the following - 
% 1. Stores and retrieves LFP for a monkey
% 2. Removes bad time segments from the data 
% 3. Computes intra probe and pairwise correlations between the two probes
% 4. Computes the distance and connectivity between two probes. 
% 5.Determines the transition channel 
% 6. Computes the mean and median pairwise correlations for different
% frequency bands for channels inside cortex. 
% 7. Plots the relationship between pairwise correlations and functional
% connectivity for different bands 
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
    chInCortexProbeAUpdated,chInCortexProbeBUpdated] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);

%% 2. Store/Retrieve the LFP (1-250 Hz) and also removing 1-5 Hz
disp(['Storing/Retrieving LFP data for ' monkeyName]);
[probe1,probe2,eeg,scalpEEG] = saveLFPDualProbe(monkeyName,hemisphere,allDates,datFileNameAll,datFileNumAll,serverPath,fs);
disp('Data Stored/Retrieved');

%% 3.Get/retrieve the distance and the connectivity values for the sites
clear distSites connSites greenMapRef 
[distSites,connVals,refSites,movSites,greenMapRef] = getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);

% Get the connectivity and distance in a vector
distSitesAll = []; connValsAll = []; 

for iDate = 1:size(allDates,1) 
    datFileNum   = datFileNumAll{iDate,1};
    distSitesAll = [distSitesAll;distSites{iDate,1}(datFileNum)];
    connValsAll  = [connValsAll; squeeze(mean(connVals{iDate,1}(:,:,datFileNum),[1,2],'omitnan'))];
end
disp(['Obtained/retrieved distance between probes and connectivity values for ' monkeyName]);

%% 4. Remove bad time segments from data and store intra probe, pairwise correlograms
 [removeTimes] = getCorrelograms(monkeyName, allDates, datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,saveFigureFlag);

%% 5. Get the transition channels for the recordings... 
% Step 1: Check powers to see if there are bad channels.
% Step 2: Eliminate bad channels and obtain the transition channels by
% computing slope of the marginals (obtained by averaging the intra-probe
% correlograms) - Mean is used here for sensitivity to outliers which is
% needed for picking the transition channel
segLen = 500; % 0.5 second segment size
winSize = 250; % overlap 
for iDate = 1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    
    for iFile = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end 
        fileNum = datFileNum(iFile);

        % Remove bad channels from Probe B
         if chInCortexProbeB{iDate}(iFile)~= 1
            probe2{fileNum,iDate}(:,badChProbeB)= [];
        end

        % Remove channels that have high impedance - from datasheets
        if strcmp(expDate,'02_21_2023') % Whiskey
            probe1{fileNum,iDate}(:,9) = [];
            probe2{fileNum,iDate}(:,8) = [];
        end

        if strcmp(expDate,'02_07_2023') % Charlie
            probe1{fileNum,iDate}(:,8) = [];
            probe2{fileNum,iDate}(:,9) = [];
        end
        
        % Remove the bad time segments
        if ~isempty(removeTimes{fileNum,iDate})
            probe1{fileNum,iDate}(removeTimes{fileNum,iDate},:) = [];
            probe2{fileNum,iDate}(removeTimes{fileNum,iDate},:) = [];
        end

        % Obtaiin data for Probe A and Probe B
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        
        if isempty(probeA) || isempty(probeB); continue; end % Check if probe A and B are not empty

        % Remove extra channels if present
        if size(probeB,2) == 33; probeB(:,1) = []; end
        if size(probeA,2) == 33; probeA(:,1) = []; end

        

        if ~strcmp(expDate,'08_08_2022');  probeB(:,badChProbeB)= []; end % Remove bad channels 14 and 22 from Probe B
      
        % Step 2: Get the intra probe marginals...
        marginalA = mean(corr(probeA,'Rows','complete'),2,'omitnan');

        fxA = movmean(gradient(marginalA),2,'omitnan'); % Get the slope of the marginals
        estChInCortexA{iDate}(fileNum,1) = find(fxA == max(fxA)); % Find the channel with maximum slope
        
        if estChInCortexA{iDate}(fileNum,1)+19>size(probeA,2) || strcmp(expDate,'10_17_2022') % Probe A on this date has spacing 50um
            estChInCortexA{iDate}(fileNum,2) = size(probeA,2);
        else
            estChInCortexA{iDate}(fileNum,2) = estChInCortexA{iDate}(fileNum,1)+19;
        end
        
        % Repeat the above process for Probe B
        marginalB = mean(corr(probeB,'Rows','complete'),2,'omitnan');
        fxB = movmean(gradient(marginalB),2,'omitnan');
        estChInCortexB{iDate}(fileNum,1) = find(fxB == max(fxB));

        if estChInCortexB{iDate}(fileNum,1)+19>size(probeB,2)
            estChInCortexB{iDate}(fileNum,2) = size(probeB,2);
        else
            estChInCortexB{iDate}(fileNum,2) = estChInCortexB{iDate}(fileNum,1)+19;
        end
    end
end 


%{ 
%Normalize data
for iDate = 1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end 
        fileNum = datFileNum(iFile);

        chA = estChInCortexA{iDate}(fileNum,:);
        chB = estChInCortexB{iDate}(fileNum,:);

        muA = mean(probe1{fileNum,iDate}(:,chA(1):chA(2)),2,'omitnan');
        stdA = std(probe1{fileNum,iDate}(:,chA(1):chA(2)),[],2,'omitnan');
        probe1{fileNum,iDate}(:,chA(1):chA(2))  = (probe1{fileNum,iDate}(:,chA(1):chA(2))- muA)./stdA; 

        muB = mean(probe2{fileNum,iDate}(:,chB(1):chB(2)),2,'omitnan');
        stdB = std(probe2{fileNum,iDate}(:,chB(1):chB(2)),[],2,'omitnan');
        probe2{fileNum,iDate}(:,chB(1):chB(2))  = (probe2{fileNum,iDate}(:,chB(1):chB(2))- muB)./stdB; 
    end
end 
%}
%% 6. Get the powers - instantaneous power, bandlimited power, spectrogram and pairwise correlations for each animal 
rowIdx = 1; 
bandLabels = {'Wideband';'Alpha band'; 'Beta band'; 'Gamma band'}; 
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
        if size(probeB,2) == 33; probeB(:,1) = []; end
        if size(probeA,2) == 33; probeA(:,1) = []; end

        % Get power - instantaneous, power spectrum and spectrogram
%         [powA{fileNum,iDate},powB{fileNum,iDate},powEEG{fileNum,iDate},specA{fileNum,iDate},specB{fileNum,iDate}...
%             ,specEEG{fileNum,iDate},instPowA{fileNum,iDate},instPowB{fileNum,iDate},freqVals,freqValsSpec,timeValsSpec] =...
%             getAllPowerData(probeA,probeB,eeg{fileNum,iDate},fs,monkeyName,expDate,fileNum);

        % Get pairwise correlations
        corrA = max(imgaussfilt(corr(probeA),1),0);
        corrB = max(imgaussfilt(corr(probeB),1),0);
        lowIntraCorr(rowIdx,1) = mean(corrA,'all')<0.29|  mean(corrB,'all')<0.29;
        
        disp('Obtaining pairwise correlations...');
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
                rowIdx = rowIdx+1;
                continue;
            end

            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate})
                meanPairCorr(rowIdx,iBand) = NaN; maxPairCorr(rowIdx,iBand) = NaN;
                medPairCorr(rowIdx,iBand)  = NaN; rowIdx = rowIdx+1;
                continue;
            end

            clear xA yA chA chB
            chA = estChInCortexA{iDate}(fileNum,:);

            if chA(1) == 0
                meanPairCorr(rowIdx,iBand)   = NaN; maxPairCorr(rowIdx,iBand)    = NaN;
                medPairCorr(rowIdx,iBand)    = NaN; meanIntraCorrA(rowIdx,iBand) = NaN;
                meanIntraCorrB(rowIdx,iBand) = NaN; medIntraCorrA(rowIdx,iBand)  = NaN;
                medIntraCorrB(rowIdx,iBand)  = NaN; rowIdx                       = rowIdx+1;
                continue;
            end

            chB = estChInCortexB{iDate}(fileNum,:);
            if chB(1) == 0
                meanPairCorr(rowIdx,iBand) = NaN;maxPairCorr(rowIdx,iBand) = NaN;
                medPairCorr(rowIdx,iBand) = NaN; rowIdx = rowIdx+1;
                continue;
            end

            switch iBand % Get the correlation vs connectivity vs distance plots for different bands...
                case 1 % Wideband
                    xA = probeA(:,chA(1):chA(2));
                    yA = probeB(:,chB(1):chB(2));

                case 2 % Alpha band
                    xA = filtfilt(bA,aA,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bA,aA,probeB(:,chB(1):chB(2)));

                case 3 % Beta band
                    xA = filtfilt(bB,aB,probeA(:,chA(1):chA(2)));
                    yA = filtfilt(bB,aB,probeB(:,chB(1):chB(2)));

                case 4 % Gamma band
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
% oneIdx = find(meanIntraCorrB(:,1) == 1);
connValsAll(unique([oneIdx;lowIntraCorr;nanIdx]))=[]; 
distSitesAll(unique([oneIdx;lowIntraCorr;nanIdx])) = [];
meanIntraCorrA(unique([oneIdx;lowIntraCorr;nanIdx]),:)=[]; 
meanIntraCorrB(unique([oneIdx;lowIntraCorr;nanIdx]),:)=[]; 
medIntraCorrA(unique([oneIdx;lowIntraCorr;nanIdx]),:)=[]; 
medIntraCorrB(unique([oneIdx;lowIntraCorr;nanIdx]),:)=[]; 

meanPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:) = []; 
medPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:)  = []; 
maxPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:)  = [];
movMeanPairCorr(unique([oneIdx;lowIntraCorr;nanIdx]),:) = [];

%%
[connSort, connIdx] = sort(connValsAll); 
meanCorrProbes_S  = meanPairCorr(connIdx,:);
meanClipA_S  = meanIntraCorrA(connIdx,:); 
meanClipB_S  = meanIntraCorrB(connIdx,:); 
medClipA_S  = medIntraCorrA(connIdx,:); 
medClipB_S  = medIntraCorrB(connIdx,:);
movMean_S   = movMeanPairCorr(connIdx,:);
medPair_S   = medPairCorr(connIdx,:); 
maxPair_S   = maxPairCorr(connIdx,:); 
dist_S     = distSitesAll(connIdx);
len        =  length(connSort);
%%
figure;
cVal = {'k';'b';'r';'m'};
for iBand = 1:4
    subplot(2,2,iBand);
    scatter(connSort,meanCorrProbes_S(:,iBand),35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand});
    xlim([-0.3 1]); ylim([-0.1 0.1]);
    xlabel('Functional Connectivity'); ylabel('Pairwise Correlations');
    hold on; box on;
    title(bandLabels{iBand})

    xFit = linspace(min(connSort),max(connSort),1000);
    pairWiseCorr = meanCorrProbes_S(:,iBand);
    pairWiseCorr(isnan(pairWiseCorr)) = 0;
    [f,gof]  = fit(connSort,pairWiseCorr,'poly1');
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    r2Final (iBand) = r2; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);

end 
%%
figure;
cVal = {'k';'b';'r';'m'};
for iBand = 1:4
    subplot(2,2,iBand);
    scatter(dist_S,meanCorrProbes_S(:,iBand),35,'filled','MarkerFaceColor',cVal{iBand},'MarkerEdgeColor',cVal{iBand});
    xlim([-0.3 1]); ylim([-0.1 0.1]);
    xlabel('Distance (mm)'); ylabel('Pairwise Correlations');
    hold on; box on;
    title(bandLabels{iBand})

    xFit = linspace(min(connSort),max(connSort),1000);
    pairWiseCorr = meanCorrProbes_S(:,iBand);
    pairWiseCorr(isnan(pairWiseCorr)) = 0;
    [f,gof]  = fit(connSort,pairWiseCorr,'poly1');
    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    r2Final (iBand) = r2; 

    plot(xFit,fitLine,'Color',cVal{iBand},'LineStyle','--','LineWidth',2);

end 


%%
cVal = 'm';
scatter(connSort,movMean_S(:,1),35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal); 
xlim([-0.3 1]); ylim([-0.3 1]); 
hold on;

xFit = linspace(min(connSort),max(connSort),1000);
pairWiseCorr = movMean_S(:,4); 
pairWiseCorr(isnan(pairWiseCorr)) = 0; 
[f,gof]  = fit(connSort,pairWiseCorr,'poly1');
c        = coeffvalues(f);
r2       = gof.rsquare;
fitLine   = c(1)*(xFit) +c(2);

plot(xFit,fitLine,'Color',cVal,'LineStyle','--','LineWidth',2);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12,'fontweight','bold'); box on;

 f = gcf; 
exportgraphics(f,'X:\Brain bag\Figures\ConnCorr_WB.eps');

%%

cVal = 'k';
scatter(dist_S,movMean_S(:,1),35,'filled','MarkerFaceColor',cVal,'MarkerEdgeColor',cVal); 
xlim([0 15]); ylim([-0.3 1]); 
hold on;
xticks(0:2:15);
xticklabels(0:2:15);

xFit = linspace(min(dist_S),max(dist_S),1000);
pairWiseCorr = connSort;%movMean_S(:,1); 
pairWiseCorr(isnan(pairWiseCorr)) = 0; 
[f,gof]  = fit(dist_S,pairWiseCorr,'poly2');
c        = coeffvalues(f);
r2       = gof.rsquare;
fitLine   = c(1)*(xFit).^2 +c(2)*xFit + c(3);

plot(xFit,fitLine,'Color',cVal,'LineStyle','--','LineWidth',2);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12,'fontweight','bold'); box on;

 f = gcf; 
exportgraphics(f,'X:\Brain bag\Figures\CorrDist.eps');



%% Plot and save correlation between instantaneous power for different bands
for iDate = 1: size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\InstPower'],'dir')
    [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\InstPower']); 
end
    for iFile = 1:length(datFileNum)
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end
        fileNum = datFileNum(iFile);
        figure;
        for iBand = 1:4
            subplot(2,2,iBand);
            imagesc(imgaussfilt(corr(instPowA{fileNum,iDate}{iBand},instPowB{fileNum,iDate}{iBand}),1)); axis image;
            colormap jet; caxis([0 1]); colorbar; 
            title([ 'Band: ' bandLabels{iBand}]);
            yticks(1:2:size(instPowA{fileNum,iDate}{iBand},2));
            xticks(1:2:size(instPowB{fileNum,iDate}{iBand},2));
        end
        
        sgtitle(strrep([expDate ' File ' num2str(fileNum)],'_','\_'));
        f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\InstPower\InstPow_' num2str(iFile) '.png'],'Resolution',300);
        close gcf;
    end
end

%% Plot the relationship between functional connectivity, pairwise correlations and distance


%% 6. Check the state of the animal under anesthesia and remove recordings if needed 
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
for iDate = 3%1:size(eeg,2)
    for iFile = 1:size(eeg,1)
        if isempty(powEEG{iFile,iDate})% || (iDate == 3) && sum(iFile == 6:8) %||((iDate == 4) || ((iDate == 5) && sum(iFile == 1:5)))
            continue; 
        else
            plot(f,10.*log10(powEEG{iFile,iDate}));  hold on; 
        end 
    end
end 

xlabel('Frequency (Hz)'); ylabel('Power (dB)'); ylim([-20 50]); 
%%
 [bG1,aG1] = butter(3,[30 55]./(fs/2),'bandpass'); % Gamma band filtering parameters
 [bG2,aG2] = butter(3,[65 90]./(fs/2),'bandpass'); % Gamma band filtering parameters 
clear alphaCharlie gamma1Charlie gamma2Charlie wideBandCharlie alphaWhiskey gamma1Whiskey gamma2Whiskey widebandWhiskey
for iM = 1:2
    clear eeg;
    switch iM
        case 1
            eeg = eegCharlie; 
        case 2
            eeg = eegWhiskey; 
    end 
    clear alphaPow gamma1Pow  gamma2Pow wideBandPow
    for iDate = 1:size(eeg,2)
        for iFile = 1:size(eeg,1)
            if isempty(eeg{iFile,iDate}) % ||((iDate == 3) && (iFile == 6 || iFile == 7))
                continue;
            else
                if (iM == 2) && (((iDate == 3) && (iFile == 6 || iFile == 7))|| ((iDate == 4) && (iFile == 1)))
                    continue; 
                end 
                if (iM == 1) && (((iDate == 4) || ((iDate == 5) && sum(iFile == 1:5))))
                    continue; 
                end 
                alphaPow(iFile,iDate)    = 10*log10(mean(pwelch(filtfilt(bA,aA,eeg{iFile,iDate}),segLen, winSize,8:12,fs),'omitnan')); 
                betaPow(iFile,iDate)     = 10*log10(mean(pwelch(filtfilt(bB,aB,eeg{iFile,iDate}),segLen, winSize,13:30,fs),'omitnan'));
                gamma1Pow(iFile,iDate)   = 10*log10(mean(pwelch(filtfilt(bG1,aG1,eeg{iFile,iDate}),segLen, winSize,30:55,fs),'omitnan')); 
                gamma2Pow(iFile,iDate)   = 10*log10(mean(pwelch(filtfilt(bG2,aG2,eeg{iFile,iDate}),segLen, winSize,65:90,fs),'omitnan'));
                wideBandPow(iFile,iDate) = 10*log10(mean(pwelch(eeg{iFile,iDate},segLen, winSize,1:120,fs),'omitnan'));
%                 instPowEEG{iFile,iDate}  = abs(hilbert(eeg{iFile,iDate})).^2;
            end
        end
    end
   switch iM
       case 1
           alphaCharlie    = alphaPow; 
           betaCharlie     = betaPow; 
           gamma1Charlie   = gamma1Pow; 
           gamma2Charlie   = gamma2Pow; 
           wideBandCharlie = wideBandPow; 
       case 2
           alphaWhiskey    = alphaPow; 
           betaWhiskey     = betaPow;
           gamma1Whiskey   = gamma1Pow; 
           gamma2Whiskey   = gamma2Pow; 
           wideBandWhiskey = wideBandPow; 
   end   
end
%
for iBand = 1:5
    clear bandCharlie bandWhiskey bandTitle
    switch iBand
        case 1
            bandCharlie = alphaCharlie;
            bandWhiskey = alphaWhiskey; 
            bandTitle   = 'Alpha band (8-12 Hz)';
        case 2
            bandCharlie = betaCharlie; 
            bandWhiskey = betaWhiskey; 
            bandTitle   = 'Beta band (13-30 Hz)'; 
        case 3
            bandCharlie = gamma1Charlie; 
            bandWhiskey = gamma1Whiskey; 
            bandTitle   = 'Gamma band (30-55 Hz)';
        case 4
            bandCharlie = gamma2Charlie;
            bandWhiskey = gamma2Whiskey;
            bandTitle   = 'Gamma band (65-90 Hz)';
        case 5
            bandCharlie = wideBandCharlie;
            bandWhiskey = wideBandWhiskey;
            bandTitle   = 'Wide band (1-120 Hz)';
    end 

    bandCharlie(bandCharlie == 0) = []; 
    bandWhiskey(bandWhiskey == 0) = [];

    figure; 
    boxplot([bandCharlie'; bandWhiskey'],[zeros(length(bandCharlie),1); ones(length(bandWhiskey),1)],'Labels',{'Charlie Sheen','Whiskey'});
    [h,p] = ttest2(bandCharlie,bandWhiskey,'VarType','unequal'); 
    disp(['P value for ' bandTitle ' is ' num2str(p)]); 
    title(['EEG Power for ' bandTitle]);
end 

%{
removeTimes = cell(size(probe1)); 
for iDate = 1:size(allDates,1) % all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];

%     if strcmp(expDate,'08_08_2022');                 continue; end
    if ~exist([saveFolder '\TimeWise'],'dir');       [~,~] = mkdir([saveFolder '\TimeWise']);      end
    if ~exist([saveFolder '\IntraProbeCorr'],'dir'); [~,~] = mkdir([saveFolder '\IntraProbeCorr']); end
    if ~exist([saveFolder '\InterProbeCorr'],'dir'); [~,~] = mkdir([saveFolder '\InterProbeCorr']); end
    if ~exist([saveFolder '\Transition'],'dir');     [~,~] = mkdir([saveFolder '\Transition']);     end

    for iFile = 1:length(datFileNum) % all recordings for a particular experiment date
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue; end 
        fileNum = datFileNum(iFile);
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if isempty(probeA) || isempty(probeB); continue; end
        if size(probeB,2) == 33; probeB(:,1) = []; end 
        if size(probeA,2) == 33; probeA(:,1) = []; end

       if ~strcmp(expDate,'08_08_2022');  probeB(:,badChProbeB)= []; end
        winLen   = 30*fs; % Bin length - 1 minute
        binWidth = 15*fs; % Sliding window length  - 30 sec
        L        = size(probeA,1);% Length of the signal collected
         
        clc;disp(['Preprocessing Ephys data for ' expDate ' File: ' num2str(fileNum)]); 
        % Get the intra probe correlation heatmap for the recordings in a
        % sliding window bin of 15 seconds and window length of 30 seconds
        ind = 1;
        clear corrDatA corrDatB
        tVals = 1:binWidth:L;

        for iW = 1:binWidth:L-winLen+1
            corrDatA(ind,:,:) = corr(probeA(iW:iW+winLen-1,:));
            corrDatB(ind,:,:) = corr(probeB(iW:iW+winLen-1,:));
            ind = ind+1;
        end

        % Correlate each time binned correlogram with the median 
        % correlogram
        sizeA    = size(probeA,2); sizeB = size(probeB,2);
        medCorrA = reshape(squeeze(median(corrDatA,1,'omitnan')),[1 sizeA*sizeA]);
        medCorrB = reshape(squeeze(median(corrDatB,1,'omitnan')),[1 sizeB*sizeB ]);

        corrDatA_R = reshape(corrDatA,[size(corrDatA,1) sizeA*sizeA ]);
        corrDatB_R = reshape(corrDatB,[size(corrDatB,1) sizeB*sizeB ]);
        
        cA = corr(corrDatA_R',medCorrA','Rows','complete');
        cB = corr(corrDatB_R',medCorrB','Rows','complete');

        % Save the histogram of the correlation values on comparing time
        % binned correlogram with the median 
        if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\TimeWise\CorrHist_' num2str(iFile) '.png'],'file') ||  saveFigureFlag
            figure; subplot(1,2,1);
            histogram(cA); xlim([ 0 1.05]); ylim([ 0 60]); xticks(0:0.1:1);yticks(0:2:60);grid on;
             xlabel('Correlations'); ylabel('Count'); title('Probe A');
            subplot(1,2,2); histogram(cB); xlim([ 0 1.05]); yticks(0:2:60);xticks(0:0.1:1);ylim([ 0 60]);grid on;
            xlabel('Correlations'); ylabel('Count'); title('Probe B');
            sgtitle(strrep([expDate ' datafile ' num2str(fileNum)],'_','\_'));
            f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\TimeWise\CorrHist_' num2str(iFile) '.png'],'Resolution',300);
            close gcf;
        end

        % Find the time segments where the correlations are less than
        % mean-2*std and remove them
        tRemove = []; timeA = []; timeB = []; uniqueTime = [];

        if any(cA<0.75)|| any(cA<(mean(cA)-5*std(cA)))
            if ~(mean(cA)-2*std(cA)<0.5)
                timeA = find((cA-(mean(cA)-2*std(cA)))<0.01);
            else
                timeA = find((cA-(mean(cA)-std(cA)))<0.01);
            end
        end

        if ~strcmp(expDate,'08_08_2022')
            if any(cB<0.75) || any(cB<(mean(cB)-5*std(cB)))
                if ~(mean(cB)-2*std(cB)<0.5)
                    timeB = find((cB-(mean(cB)-2*std(cB)))<0.01);
                else
                    timeB = find((cB-(mean(cB)-std(cB)))<0.01);
                end
            end
        end
        
        if ~(isempty(timeA) && isempty(timeB))
            uniqueTime = tVals(unique([timeA ;timeB]));

            for iL = 1:length(uniqueTime)
                tRemove = [tRemove uniqueTime(iL):uniqueTime(iL)+winLen-1];
            end
            removeTimes{fileNum,iDate} = tRemove;
        end

        % Save one minute bin (60 second sliding window) intra probe
        % correlations for reference
        ind = 1; 
        clear corrDatAOneMin corrDatBOneMin; 
        for iW = 1:60*fs:L-(60*fs)+1
            corrDatAOneMin(ind,:,:) = corr(probeA(iW:iW+60*fs-1,:));
            corrDatBOneMin(ind,:,:) = corr(probeB(iW:iW+60*fs-1,:));
            ind = ind+1;
        end

        for iProbe = 1:2
            clear probeName plotDat
            switch iProbe
                case 1
                    probeName = 'Probe A';
                    plotDat   = corrDatAOneMin;
                case 2
                    probeName = 'Probe B';
                    plotDat   = corrDatBOneMin;
            end

            tValsOneMin = 0:60*fs:L;

            if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\TimeWise\IntraProbe' probeName(end) '_' num2str(iFile) '.png'],'file') ||  saveFigureFlag
                figure('units','normalized','outerposition',[0 0 1 1]);
                for iL = 1:size(plotDat,1)
                    subplot(ceil(size(plotDat,1)/5),5,iL);
                    imagesc(imgaussfilt(squeeze(plotDat(iL,:,:)),1)); colormap jet; title([num2str(tValsOneMin(iL)/1e3) ' to ' num2str(tValsOneMin(iL+1)/1e3) ' s']);
                    caxis([0 1]); if iL == 1; colorbar; end
                end
                sgtitle(strrep([expDate ' datafile ' num2str(datFileNum(iFile)) ' ' probeName ],'_','\_'));
                f = gcf;
                exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\TimeWise\IntraProbe' probeName(end) '_' num2str(iFile) '.png'],'Resolution',300);
                close gcf;
            end
        end

        % Get the intra probe correlograms (verticals and heatmaps) after
        % removing the bad time segments 
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        if size(probeA,2) == 33; probeA(:,1) = []; end 
        if size(probeB,2) == 33; probeB(:,1) = []; end 

        % Remove the bad time segments 
        if ~isempty(removeTimes{fileNum,iDate})
            probeA(removeTimes{fileNum,iDate},:) = [];
            probeB(removeTimes{fileNum,iDate},:) = [];
        end

        chA = chInCortexProbeA{iDate}(iFile)+10; % 1mm from the transition point listed in experiment notes
        if ~strcmp(expDate,'08_08_2022');chB = chInCortexProbeB{iDate}(iFile)+10; end

        for iBand = 1:3
            clear pA pB figTitle b a
            switch iBand 
                case 1
                    pA = probeA;
                    pB = probeB;
                    figTitle = 'WB';

                case 2
                    pA = filtfilt(bG,aG,probeA);
                    pB = filtfilt(bG,aG,probeB);
                    figTitle = 'Gamma band';

                case 3
                    pA = filtfilt(bA,aA,probeA);
                    pB = filtfilt(bA,aA,probeB);
                    figTitle = 'Alpha band';
            end 

            clear cortexChA chOutA chDeepA cortexChB chOutB chDeepB 
            cortexChA = pA(:,chA:chA+2);               % Correlating with channels that we know is inside the cortex
            chOutA    = pA(:,chOutCortex);             % Correlating with channels that we know is out of cortex
            chDeepA   = pA(:,chDeep);                  % Correlating with the deepest channels

            if ~strcmp(expDate,'08_08_2022') % Condition where there is only one channel for a site (applicable for Whiskey)
                cortexChB = pB(:,chB:chB+2);
                chOutB    = pB(:,chOutCortex);
                chDeepB   = pB(:,chDeep);
            else
                cortexChB = pB;
                chOutB    = zeros(size(pB));
                chDeepB   = pB;
            end

            if ~strcmp(expDate,'08_08_2022'); pB(:,badChProbeB)= []; end % Remove bad channels from probe B

            % Plot and save the intra probe correlation heatmaps...
            if ~exist([saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'file') || ~exist([saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.eps'],'file') ||  saveFigureFlag
                figure;
                imagesc(imgaussfilt(corr(pA,'Rows','complete'),1)); colormap jet; xticks(1:32); yticks(1:32); caxis([ 0 1]); colorbar;
                title(strrep([expDate ' Datafile ' num2str(iFile) ' Probe A - ' figTitle ],'_','\_'));
                f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.eps'],'ContentType','vector');
                close gcf;
            end            

            if ~exist([saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'file')|| ~exist([saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.eps'],'file') ||  saveFigureFlag
                figure;
                imagesc(imgaussfilt(corr(pB,'Rows','complete'),1)); colormap jet; xticks(1:30); yticks(1:30); caxis([ 0 1]); colorbar;
                xticklabels(probeBList); yticklabels(probeBList);title(strrep([expDate ' Datafile ' num2str(iFile) ' Probe B - ' figTitle  ],'_','\_'));
                f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.eps'],'ContentType','vector');
                close gcf;
            end

            % Get the intra probe marginals... 
            marginalA = mean(corr(pA,'Rows','complete'),2);
            marginalB = mean(corr(pB,'Rows','complete'),2);
            
            if iBand == 1
                fxA = movmean(gradient(marginalA),2);
                fxB = movmean(gradient(marginalB),2);
                estChInCortexA{iDate}(fileNum,1) = find(fxA == max(fxA));
               
                if estChInCortexA{iDate}(fileNum,1)+19>size(probeA,2)
                    estChInCortexA{iDate}(fileNum,2) = size(probeA,2);
                else
                    estChInCortexA{iDate}(fileNum,2) = estChInCortexA{iDate}(fileNum,1)+19;
                end

                estChInCortexB{iDate}(fileNum,1) = find(fxB == max(fxB));

                if estChInCortexB{iDate}(fileNum,1)+19>size(probeB,2)
                    estChInCortexB{iDate}(fileNum,2) = size(probeB,2);
                else
                    estChInCortexB{iDate}(fileNum,2) = estChInCortexB{iDate}(fileNum,1)+19;
                end
            end 
            

            if ~exist([saveFolder '\IntraProbeCorr\Marginal_File_' num2str(fileNum) '_' figTitle '.png'],'file') || ~exist([saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.eps'],'file') ||  saveFigureFlag
                figure; subplot(1,2,1); 
                
                imagesc(imgaussfilt(marginalA,1)); colormap jet; yticks(1:32); caxis([ 0 max(marginalA)]); colorbar; 
                title('Probe A');

                subplot(1,2,2); 
                imagesc(imgaussfilt(marginalB,1)); colormap jet; yticks(1:30);yticklabels(probeBList); caxis([ 0 max(marginalB)]); colorbar; 
                title('Probe B');
                sgtitle(strrep([expDate ' Datafile ' num2str(iFile) ' ' figTitle ],'_','\_'));

                f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\Marginal_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                exportgraphics(f,[saveFolder '\IntraProbeCorr\Marginal_File_' num2str(fileNum) '_' figTitle '.eps'],'ContentType','vector');
                close gcf;
            end            

            % Plot and save the intra probe verticals...
            if ~exist([saveFolder '\Transition\Transition_' num2str(fileNum) '_Probe A_' figTitle '.png'],'file') || ~exist([saveFolder '\Transition\Transition_' num2str(fileNum) '_Probe A_' figTitle '.eps'],'file') ||  saveFigureFlag
                if iBand~=3
                    for iProbe = 1:2
                        clear corrIn corrOut corrDeep probeTitle in
                        switch iProbe
                            case 1
                                corrIn = corr(pA,cortexChA);
                                corrOut = corr(pA,chOutA) ;
                                corrDeep = corr(pA,chDeepA);
                                in       = chA:chA+2;
                                probeTitle = 'Probe A';
                            case 2
                                corrIn = corr(pB,cortexChB);
                                corrOut = corr(pB,chOutB) ;
                                corrDeep = corr(pB,chDeepB);
                                in       = chB:chB+2;
                                probeTitle = 'Probe B';
                        end

                        figure;

                        subplot(1,3,1); imagesc(movmean(mean(corrIn,2,'omitnan'),5));  caxis([ 0 1]); colormap jet;
                        ylabel('Mean correlation'); title(['Ch: ' num2str(in)]);  xticks([]); yticks(1:size(corrIn,1)); colorbar;
                        if iProbe == 2; yticks(1:30); yticklabels(probeBList);end

                        subplot(1,3,2); imagesc(movmean(mean(corrOut,2,'omitnan'),5)); caxis([ 0 1]); colormap jet;
                        ylabel('Mean correlation');title( ['Ch: ' num2str(chOutCortex)]);xticks([]); yticks(1:size(corrOut,1)); colorbar;
                        if iProbe == 2; yticks(1:30); yticklabels(probeBList);end

                        subplot(1,3,3); imagesc(movmean(mean(corrDeep,2,'omitnan'),5));    caxis([0 1]); colormap jet;
                        ylabel('Mean correlation'); title(['Ch: ' num2str(chDeep)]);xticks([]); yticks(1:size(corrDeep,1)); colorbar;
                        if iProbe == 2; yticks(1:30); yticklabels(probeBList);end

                        sgtitle(strrep([expDate ' Datafile ' num2str(fileNum) ' ' probeTitle ' - ' figTitle  ],'_','\_'));

                        f = gcf; exportgraphics(f,[saveFolder '\Transition\Transition_' num2str(fileNum) '_' probeTitle '_' figTitle '.png'],'Resolution',300);
                        exportgraphics(f,[saveFolder '\Transition\Transition_' num2str(fileNum) '_' probeTitle '_' figTitle '.eps'],'ContentType','vector');
                        close gcf;

                    end
                end
            end
            
            % Plot and save the pairwise correlation heatmaps
            if ~exist([saveFolder '\InterProbeCorr\AllCh_File_' num2str(fileNum) '_' figTitle '.png'],'file') || ~exist([saveFolder '\InterProbeCorr\AllCh_File_' num2str(fileNum) '_' figTitle '.eps'],'file') ||  saveFigureFlag
                figure; imagesc(imgaussfilt(corr(pA,pB,'rows','complete'),1)); colormap jet; caxis([ 0 1]);
                title(strrep([expDate ' Datafile ' num2str(fileNum) ': Pairwise correlation -  ' figTitle],'_','\_')); xlabel('Probe B'); ylabel('Probe A'); colorbar;
                yticks(1:32); xticks(1:30); xticklabels(probeBList);
                f = gcf; exportgraphics(f,[saveFolder '\InterProbeCorr\AllCh_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                exportgraphics(f,[saveFolder '\InterProbeCorr\AllCh_File_' num2str(fileNum) '_' figTitle '.eps'],'ContentType','vector');
                close gcf;
            end
        end
    end
end
clc;
disp(['Preprocessing complete for ' monkeyName]);
%}
%% Determine recordings that are bad and are to be rejected from further analysis
% Find the mean intra probe correlations for all 32 channels and remove the
% recordings where either Probe A or Probe B has intra probe correlations
% less than 0.29
clear corrMaxProbeA corrMaxProbeB corrMeanProbeA corrMeanProbeB meanClipA meanClipB medClipA medClipB
if ~exist('plotCorrFigs','var') || isempty(plotCorrFigs); plotCorrFigs = 1;    end 
rowIdx = 1; 

for iDate = 1:size(allDates,1)
    clear expDate datFileNum
      expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum) % Recording #
        clear corrChA corrChAOut corrChARandom corrChB corrChBOut corrChBRandom maxA maxB
        fileNum = datFileNum(iFile);
        if strcmp(expDate,'09_19_2022') && fileNum == 4; rowIdx = rowIdx+1; continue; end 
         if strcmp(expDate,'10_17_2022') && (fileNum == 6 || fileNum == 7 || fileNum == 8); rowIdx = rowIdx+1; continue; end  

        % Get the LFP of the probes 
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        if size(probeB,2)>32 % Finding channels that have colums of zero
            idx = find(sum(probeB == 0));
            if ~isempty(idx)
                probeB(:,idx) = [];
            end
        end

        if isempty(probeA) && isempty(probeB); rowIdx = rowIdx+1;  continue; end % To continue if there is no LFP data
%         if ~strcmp(expDate,'08_08_2022'); probeB(:,badChProbeB)= []; end % Remove bad channels from probe B

        % Remove the bad time segments
        if ~isempty(removeTimes{fileNum,iDate})
            probeA(removeTimes{fileNum,iDate},:) = [];
            probeB(removeTimes{fileNum,iDate},:) = [];
        end

         corrA               = imgaussfilt(corr(probeA),1); corrB = imgaussfilt(corr(probeB),1); 
         meanA(rowIdx,1)     = mean(corrA,'all');          meanB(rowIdx,1) = mean(corrB,'all');
         medA(rowIdx,1)      = median(corrA,'all');        medB(rowIdx,1) = median(corrB,'all');
         corrA               = max(corrA,0);               corrB = max(corrB,0);
         meanClipA(rowIdx,1) = mean(corrA,'all');         meanClipB(rowIdx,1) = mean(corrB,'all');
         medClipA(rowIdx,1)  = median(corrA,'all');       medClipB(rowIdx,1) = median(corrB,'all');
         rowIdx = rowIdx+1; 
    end
end
lowIntraCorr =[];% find(meanClipA<0.29| meanClipB<0.29);  % Threshold for Charlie/Whiskey

%{
% Begin determining transition for all recordings... 
fileIdx = 1; 
clear estChInCortexA estChInCortexB
for iDate = 1:size(allDates,1) % Experiment date
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 1:length(datFileNum) % Recording #
        clear corrChA corrChAOut corrChARandom corrChB corrChBOut corrChBRandom maxA maxB
        fileNum = datFileNum(iFile);

        if any(lowIntraCorr == fileIdx)
            fileIdx = fileIdx+1; 
            estChInCortexA{iDate}(fileNum,:)= [0 0];
            estChInCortexB{iDate}(fileNum,:) = [0 0];
            continue;
        end

%         if strcmp(expDate,'09_19_2022') && fileNum == 4; fileIdx = fileIdx+1; continue; end 
        if strcmp(expDate,'10_17_2022') && (fileNum == 6 || fileNum == 7 || fileNum == 8); fileIdx = fileIdx+1; continue; end  

        % Get the LFP of the probes 
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if size(probeB,2)>32 % Finding channels that have colums of zero
            idx = find(sum(probeB == 0));
            if ~isempty(idx)
                probeB(:,idx) = [];
            end
        end

        if isempty(probeA) && isempty(probeB); fileIdx = fileIdx+1; continue; end % To continue if there is no LFP data

         % Remove the bad time segments
        if ~isempty(removeTimes{fileNum,iDate})
            probeA(removeTimes{fileNum,iDate},:) = [];
            probeB(removeTimes{fileNum,iDate},:) = [];
        end

        chA       = chInCortexProbeA{iDate}(iFile)+10; % 1mm from the transition point listed in experiment notes
        cortexChA = probeA(:,chA:chA+2);               % Correlating with channels that we know is inside the cortex
        chOutA    = probeA(:,chOutCortex);             % Correlating with channels that we know is out of cortex
        chDeepA   = probeA(:,chDeep);                  % Correlating with the deepest channels

        if ~strcmp(expDate,'08_08_2022') % Condition where there is only one channel for a site (applicable for Whiskey)
            chB       = chInCortexProbeB{iDate}(iFile)+10;
            cortexChB = probeB(:,chB:chB+2);
            chOutB    = probeB(:,chOutCortex);
            chDeepB   = probeB(:,chDeep);
        else
            cortexChB = probeB;
            chOutB    = zeros(size(probeB));
            chDeepB   = probeB;
        end

        if ~strcmp(expDate,'08_08_2022'); probeB(:,badChProbeB)= []; end % Remove bad channels from probe B

        % Condition 1: To check if the probe is too deep. Check if atleast
        % the first 8 contacts have correlation value greater than 0.65 and
        % or the median to be greater than 0.7. If yes, Ch 1 is the
        % superficial channel and channel 20 is the deepest.
        probeADeepFlag = 0;
        probeBDeepFlag = 0;
        
        % Probe A 
        corrA  = mean(corr(probeA,chOutA),2,'omitnan');
        if median(corrA)>0.7 || (sum(corrA(1:10)>0.65)>8) 
            estChInCortexA{iDate}(fileNum,1) = 1; estChInCortexA{iDate}(fileNum,2) = 20;
            probeADeepFlag  = 1;
        end

        % Probe B
        if ~strcmp(expDate,'08_08_2022')
            corrB    = mean(corr(probeB,chOutB),2,'omitnan');
            if median(corrB)>0.7 ||(sum(corrB(1:10)>0.65)>8 )
                estChInCortexB{iDate}(fileNum,1) = 1; estChInCortexB{iDate}(fileNum,2) = 20;
                probeBDeepFlag  = 1;
            end
        else
            estChInCortexB{iDate}(fileNum,1) = 1; estChInCortexB{iDate}(fileNum,2) = 1;
            probeBDeepFlag = 1;
        end

        % Condition 2: If the probe is not too deep, get the slope of the
        % correlations for channels 1mm deep and also for the deepest
        % channels. Find the channel that contains the maximum slope. If
        % If the correlation of the transition channel is less than the
        % median-std, then the process is repeated till the optimal
        % transition is found. Once found, the average of the two
        % transition channels are taken as the transition point for the
        % probe

        % Probe A
        if ~probeADeepFlag
            corrA =mean(corr(probeA,cortexChA),2,'omitnan'); stdA  = std(corrA); medInA = median(corrA);
            fxA    = movmean(gradient(corrA),2); % Transition determined from wide band
            maxAcortex = find(fxA == max(fxA(2:end)));

            corrADeep = mean(corr(probeA,chDeepA),2,'omitnan');
            fxADeep    = movmean(gradient(corrADeep),2); % Transition determined from wide band
            maxAdeep = find(fxADeep == max(fxADeep(2:end)));
            avgA = round(mean([maxAcortex,maxAdeep],'omitnan'));

            while ~((abs(corrA(avgA)-(medInA-stdA)))<0.01) && (corrA(avgA)<=(medInA-stdA))
                maxAcortex = find(fxA == max(fxA(maxAcortex+1:end)));
                maxAdeep = find(fxADeep == max(fxADeep(maxAdeep+1:end)));
                avgA = round(mean([maxAcortex,maxAdeep],'omitnan'));
            end 

            estChInCortexA{iDate}(fileNum,1) = avgA;
            if avgA+19>size(probeA,2)
                estChInCortexA{iDate}(fileNum,2) = size(probeA,2);
            else
                estChInCortexA{iDate}(fileNum,2) = avgA+19;
            end
        end

        % Probe B
        if ~probeBDeepFlag
            corrB = mean(corr(probeB,cortexChB),2,'omitnan'); stdB = std(corrB); medInB = median(corrB);
            fxB   = movmean(gradient(corrB),2); % Transition determined from wide band
            maxBcortex = find(fxB == max(fxB(2:end)));

            corrBDeep = mean(corr(probeB,chDeepB),2,'omitnan');
            fxBDeep    = movmean(gradient(corrBDeep),2); % Transition determined from wide band
            maxBdeep = find(fxBDeep == max(fxBDeep(2:end)));
            avgB = round(mean([maxBcortex,maxBdeep],'omitnan'));

            while ~((abs(corrB(avgB)-(medInB-stdB)))<0.01) && (corrB(avgB)<=(medInB-stdB))
                maxBcortex = find(fxB == max(fxB(maxBcortex+1:end)));
                maxBdeep = find(fxBDeep == max(fxBDeep(maxBdeep+1:end)));
                avgB = round(mean([maxBcortex,maxBdeep],'omitnan'));
            end
           
            estChInCortexB{iDate}(fileNum,1) = avgB;
            if avgB+19>size(probeB,2)
                estChInCortexB{iDate}(fileNum,2) = size(probeB,2);
            else
                estChInCortexB{iDate}(fileNum,2) = avgB+19;
            end
        end   
        fileIdx = fileIdx+1; 
    end
end
%}
%% 6. Get the correlation vs connectivity vs distance for WB,alpha and gamma frequencies
saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results'];
clear nanIdx negIdxCorr negIdxConn  meanPairCorr medPairCorr maxPairCorr
if ~exist([saveFolder '\Corr_Conn_Dist'],'dir'); [~,~] = mkdir([saveFolder '\Corr_Conn_Dist']); end

% Get intra probe correlations and pairwise correlations after taking the channels in the cortex...  
for iBand = 1:4
    disp(['Computing measures of central tendency for the band: ' bandLabels{iBand}]); 
    rowIdx = 1;
    for iDate = 1:size(allDates,1)
        clear expDate datFileNum
        expDate = allDates(iDate,:);
        datFileNum = datFileNumAll{iDate,1};
        for iFile = 1:length(datFileNum)
            clear fileNum
            fileNum = datFileNum(iFile);
            if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7))||(strcmp(expDate,'02_07_2023') && (fileNum == 2 )) ||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))||(strcmp(expDate,'09_19_2022') && (fileNum == 4))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) || 
                meanPairCorr(rowIdx,iBand) = NaN; 
                maxPairCorr(rowIdx,iBand)  = NaN; 
                medPairCorr(rowIdx,iBand)  = NaN;

                meanIntraCorrA(rowIdx,iBand) = NaN; 
                meanIntraCorrB(rowIdx,iBand) = NaN; 
                medIntraCorrA(rowIdx,iBand)  = NaN; 
                medIntraCorrB(rowIdx,iBand)  = NaN; 
                rowIdx = rowIdx+1;
                continue;
            end
            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate}); meanPairCorr(rowIdx,iBand) = NaN;maxPairCorr(rowIdx,iBand) = NaN; 
                medPairCorr(rowIdx,iBand) = NaN; rowIdx = rowIdx+1; continue; 
            end

            clear xA yA chA chB
            chA(1) = estChInCortexA{iDate}(fileNum);
%             chA(1) = chInCortexProbeA{iDate}(iFile);
            if chA(1) == 0
                meanPairCorr(rowIdx,iBand) = NaN; maxPairCorr(rowIdx,iBand) = NaN;
                medPairCorr(rowIdx,iBand) = NaN;  meanIntraCorrA(rowIdx,iBand) = NaN;
                meanIntraCorrB(rowIdx,iBand) = NaN; medIntraCorrA(rowIdx,iBand) = NaN;
                medIntraCorrB(rowIdx,iBand) = NaN; rowIdx = rowIdx+1;
                continue;
            end
            if chA(1)+19 > 32
                chA(2) = 32; 
            else
                chA(2) = chA(1)+19; 
            end
           
            xA = probe1{fileNum,iDate}(:,chA(1):chA(2));

            if chInCortexProbeB{iDate}(iFile)~= 1
                yA = probe2{fileNum,iDate};  if strcmp(expDate,'11_01_2021'); yA = yA(:,2:end); end 
                yA(:,badChProbeB)= [];
                chB(1) =estChInCortexB{iDate}(fileNum);
%                 chB(1) = chInCortexProbeB{iDate}(iFile);

                if chB(1) == 0; meanPairCorr(rowIdx,iBand) = NaN;maxPairCorr(rowIdx,iBand) = NaN;
                    medPairCorr(rowIdx,iBand) = NaN; rowIdx = rowIdx+1; continue;  end
                if chB(1)+19 > size(yA,2)
                    chB(2) = size(yA,2);
                else
                    chB(2) = chB(1)+19;
                end
                yA = yA(:,chB(1):chB(2));
            else
                yA = probe2{fileNum,iDate};
                if strcmp(expDate,'11_01_2021'); yA = yA(:,2:end); end
            end
            if  sum(isnan(yA),'all') >= size(yA,1) || sum(isnan(xA),'all') >= size(xA,1); meanPairCorr(rowIdx,iBand) = NaN; maxPairCorr(rowIdx,iBand) = NaN;
                    medPairCorr(rowIdx,iBand) = NaN; rowIdx = rowIdx+1; continue; end 
            switch iBand % Get the correlation vs connectivity vs distance plots for different bands...
                case 1 % WB 
                    xA = xA; yA = yA; 

                case 2 % Gamma band
                    xA = filtfilt(bG,aG,xA);
                    yA = filtfilt(bG,aG,yA);

                case 3 % Alpha band
                    xA = filtfilt(bA,aA,xA);
                    yA = filtfilt(bA,aA,yA);

                case 4 % Beta band
                    xA = filtfilt(bB,aB,xA); 
                    yA = filtfilt(bB,aB,yA);
            end

            % Remove the bad time segments... 
            if ~isempty(removeTimes{fileNum,iDate})
                xA(removeTimes{fileNum,iDate},:) = [];
                yA(removeTimes{fileNum,iDate},:) = [];
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

%             [fA,xAFreq] = getFFT(xA,fs,0); 
%             [fB,yAFreq] = getFFT(yA,fs,0);
            % Get moving window pairwise correlations between the two
            % probes...
            ind = 1; 
            L = size(xA,1); 
            clear movMeanPairCorrTemp
            for iW = 1:30*fs:L-(60*fs)+1
                movMeanPairCorrTemp(ind,1) = mean(corr(xA(iW:iW+60*fs-1,:),yA(iW:iW+60*fs-1,:)),'all','omitnan');
                ind = ind+1;
            end
            movMeanPairCorr(rowIdx,iBand) = mean(movMeanPairCorrTemp,1,'omitnan');
             rowIdx = rowIdx+1;
        end
    end
end
% lowIntraCorr = find(meanIntraCorrA(:,1)<0.29| meanIntraCorrB(:,1)<0.29);
%%
figure; 
colororder({'k','r'});
yyaxis left; 
plot(1:len,connSort,'k.','MarkerSize',12);  
ylim([-1 1]); 
yticks(-1:0.2:1); 
ylabel('Connectivity'); 
xlim([0 len+1]); 

yyaxis right; 
plot(1:len,movMean_S(:,2),'r*','MarkerSize',12); 
xlabel('Count'); ylabel('Pairwise correlation'); ylim([-1 1]); yticks(-1:0.2:1); 
 hold on; %grid on;

xFit     = linspace(1,len,1000);
pairWiseCorr = movMean_S(:,1); 
pairWiseCorr(isnan(pairWiseCorr)) = 0; 
[f,gof]  = fit((1:len)',pairWiseCorr,'poly1');
c        = coeffvalues(f);
r2       = gof.rsquare;
fitLine   = c(1)*(xFit) +c(2);

plot(xFit,fitLine,'Color','r','LineStyle','--','LineWidth',1);

xticks(0:len);
%%
yyaxis right; 
hold on; 
plot(1:len, meanClipA_S(:,1),'ro','MarkerSize',5); 
plot(1:len, meanClipB_S(:,1),'ko','MarkerSize',5); 
ylabel('Correlation values');
legend('Connectivity','Pairwise correlation','Fit for pairwise correlation','Probe A','Probe B','location','southwest');

%%
figure; 
plot(meanPairCorr(:,1),meanIntraCorrA(:,1),'ro','MarkerSize',5); 
hold on; 
plot(meanPairCorr(:,1),meanIntraCorrB(:,1),'bo','MarkerSize',5); 
xlabel('Inter probe correlations');ylabel('Intra probe correlations');
grid on; xlim([ -0.5 1]); ylim([ -0.5 1]); xticks(-0.5:0.1:1); yticks(-0.5:0.1:1);
legend('Probe A','Probe B','location','southwest');

%% Same plot, but with fit for intra probes... 
figure; 
colororder({'k','k'});
yyaxis left; 
plot(1:len,meanCorrProbes_S(:,1),'b*','MarkerSize',12); 
xlabel('Count'); ylabel('Pairwise correlation'); ylim([-1 1]); yticks(-1:0.1:1); 
xticks(1:len);xlim([0 len+1]);
grid on; hold on; 

xFit     = linspace(1,len,1000);
pairWiseCorr = meanCorrProbes_S(:,1); 
pairWiseCorr(isnan(pairWiseCorr)) = 0; 
[f,gof]  = fit((1:len)',pairWiseCorr,'poly1');
c        = coeffvalues(f);
r2       = gof.rsquare;
fitLine   = c(1)*(xFit) +c(2);

plot(xFit,fitLine,'Color','b','LineStyle','--','LineWidth',1);

yyaxis right; 
plot(1:len,meanClipA_S(:,1),'ro','MarkerSize',5); 
hold on; ylim([-1 1]); yticks(-1:0.1:1);
xFit     = linspace(1,len,1000);
pairWiseCorr = meanClipA_S(:,1); 
pairWiseCorr(isnan(pairWiseCorr)) = 0; 
[f,gof]  = fit((1:len)',pairWiseCorr,'poly1');
c        = coeffvalues(f);
r2       = gof.rsquare;
fitLine   = c(1)*(xFit) +c(2);

plot(xFit,fitLine,'Color','r','LineStyle','--','LineWidth',1);

plot(1:len,meanClipB_S(:,1),'ko','MarkerSize',5); 
xFit     = linspace(1,len,1000);
pairWiseCorr = meanClipB_S(:,1); 
pairWiseCorr(isnan(pairWiseCorr)) = 0; 
[f,gof]  = fit((1:len)',pairWiseCorr,'poly1');
c        = coeffvalues(f);
r2       = gof.rsquare;
fitLine   = c(1)*(xFit) +c(2);

plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',1);

ylabel('Intra probe correlations');
grid on; hold on; 
legend('Pairwise','Pairwise- fit','Probe A','Probe A-fit','Probe B','Probe B-fit','location','southwest')

xticks(0:len);

%%
nanIdx = find(isnan(meanPairCorr(:,1))); 
meanPairCorr(nanIdx,:) = []; 
medPairCorr(nanIdx,:)  = [];
maxPairCorr(nanIdx,:)  = [];
connValsAll(nanIdx)      = []; 
distSitesAll(nanIdx)     = []; 

siteConnAllTemp = cell2mat(siteConnAll)'; 
siteConnAllTemp(nanIdx) = [];
siteConnAllTemp = logical(siteConnAllTemp);

%%  Determine the transition channel - to modify this... 
% [estChInCortexA,estChInCortexB] = getTransitionChannels(monkeyName,hemisphere,allDates,datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,1);


%% Plot the correlation vs connectivity vs distance for different frequency bands 
distSitesTemp      = distSitesAll; 
connValsTemp       = connValsAll;
meanCorrProbesTemp = meanPairCorr; 
for iClip = 1:2
    switch iClip
        case 1
            clipName = 'noclip';
        case 2
            clipName = 'clip';
    end
    for iBand = 1:3
        switch iBand
            case 1
                bandName = 'WB';
            case 2
                bandName = 'Gamma';
            case 3
                bandName = 'Alpha';
        end
        if iClip == 2
            % Remove the negative correlations and connectivity values
            negIdxCorr = meanCorrProbesTemp(:,iBand)<0;
            negIdxConn = connValsTemp<0;

            meanCorrProbesTemp(negIdxCorr,iBand) = 0;
            connValsTemp(negIdxConn)    = 0;
            yaxesLim = [0 1];
        else
            yaxesLim = [-1 1];
        end

        % Plotting 2D plots of correlation vs connectivity
        figure;
        plot(connValsTemp,meanCorrProbesTemp(:,iBand),'k.','MarkerSize',10); hold on; grid on;
        xlabel('Connectivity'); ylabel('Correlation'); title(['Mean Correlation vs connectivity - ' bandName ]);
        xlim([min(connValsTemp)-0.1 max(connValsTemp)+0.1]); ylim(yaxesLim);
        xticks(round(min(connValsTemp)-0.1:0.1:max(connValsTemp)+0.1,1)); yticks(yaxesLim(1):0.1:yaxesLim(2));

%         xFit     = linspace(min(connValsTemp),max(connValsTemp),1000);
%         [f,gof]  = fit(connValsTemp,meanCorrProbesTemp(:,iBand),'poly1');
%         c        = coeffvalues(f);
%         r2       = gof.rsquare;
%         fitLine   = c(1)*(xFit) +c(2);
% 
%         plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',1);
%         text(max(connValsTemp),0.7,['R^2: ' num2str(r2)],'Color','k');
        f = gcf; exportgraphics(f,[saveFolder '\Corr_Conn_Dist\MeanCorr_Conn_' clipName '_' bandName '.png'],'Resolution',300); close gcf;

        % Plotting 2D plots of correlation vs distance
        figure;
        for iLine = 1:2
            clear allVals corrVals colorVals
            switch iLine
                case 1
                    allVals   = [distSitesTemp(siteConnAllTemp)];
                    corrVals  = [ meanCorrProbesTemp(siteConnAllTemp,iBand)];
                    colorVals = 'r';
                case 2
                    allVals = [distSitesTemp(~siteConnAllTemp)];
                    corrVals =[ meanCorrProbesTemp(~siteConnAllTemp,iBand)];
                    colorVals = 'b';
            end

            plot(allVals,corrVals,'.','Color',colorVals,'MarkerSize',12);grid on; hold on;
            xlabel('Distance (mm)'); ylabel('Correlation'); title(['Mean Correlation vs distance - ' bandName ]);
        end
        xlim([0 max(distSitesTemp)+1]); ylim(yaxesLim); 
        xticks(1:max(distSitesTemp)+1); yticks(yaxesLim(1):0.1:yaxesLim(2));
        legend('Connected sites','Not connected sites','Location','northeast');
        f = gcf; exportgraphics(f,[saveFolder '\Corr_Conn_Dist\MeanCorr_Dist_' clipName '_' bandName '.png'],'Resolution',300); close gcf;
    end
end

%% Correlation vs connectivity dependent on distance
% Take the correlation and connectivity for certain distance values:0-3;
% 3-6; 6-9; 9-12; 12-15 mm
saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\Results'];
if ~exist([saveFolder '\Corr_Conn_Dist\DistConstrain'],'dir'); [~,~] = mkdir([saveFolder '\Corr_Conn_Dist\DistConstrain']); end
distSitesTemp      =  distSitesAll;%distUpdated;%distSitesAll;
connValsTemp       = connValsAll;%connUpdated;%connValsAll;
meanCorrProbesTemp = meanPairCorr;% meanCorrUpdated;%maxCorrProbes;
maxVal = round(max(distSitesTemp));
for iClip = 1:2
    switch iClip
        case 1
            clipLabel ='no clip';
            saveLabel ='noclip';
        case 2
            clipLabel = 'Clipped at 0';
            saveLabel = 'clip';
    end
    for iDist = 0:3:maxVal
        findDist = find(distSitesTemp>iDist & distSitesTemp<=iDist+3);
        tempDist = distSitesTemp(findDist);
        tempConn = connValsTemp(findDist);
        tempMeanCorr = meanCorrProbesTemp(findDist,:);
        figure;
        for iBand = 1%:3
            switch iBand
                case 1
                    bandCol = 'k';
                    bandName = 'WB';
                case 2
                    bandCol = 'r';
                    bandName ='Gamma';
                case 3
                    bandCol = 'b';
                    bandName = 'Alpha';
            end

            if iClip == 2
                % Remove the negative correlations and connectivity values
                negIdxCorr = tempMeanCorr(:,iBand)<0;
                negIdxConn = tempConn<0;

                tempMeanCorr(negIdxCorr,iBand) = 0;
                tempConn(negIdxConn)    = 0;
                locVal  = 'northwest';
            else
                locVal = 'southeast';
            end

            plot(tempConn,tempMeanCorr(:,iBand),'.','Color',bandCol,'MarkerSize',15, ...
                'DisplayName',bandName); hold on; grid on;
            xlabel('Connectivity'); ylabel('Correlation'); title(['Mean Correlation vs connectivity - ' clipLabel ' - Distance: ' num2str(iDist) ' to ' num2str(iDist+3) ' mm']);
            xlim([min(tempConn)-0.1 max(tempConn)+0.1]); ylim([-1 1]);

            if length(tempConn)>2
                xFit = linspace(min(tempConn),max(tempConn),1000);
                [f,gof]  = fit(tempConn,tempMeanCorr(:,iBand),'poly1');
                c        = coeffvalues(f);
                r2       = gof.rsquare;
                fitLine   = c(1)*(xFit) +c(2);

                plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',1);
                text(max(tempConn) - 0.1, 0.9,['R^2 : ' num2str(r2)]);
            end
            
        end
        legend('Location',locVal);
        f = gcf; exportgraphics(f,[saveFolder '\Corr_Conn_Dist\DistConstrain\WB_MeanCorr_Conn_' saveLabel '_Dist_ ' num2str(iDist) '_' num2str(iDist+3) '.png'],'Resolution',300); close gcf;
    end
end


%% 6.  Get the correlation heatmaps and pairwise correlations for the recordings
clear corrAll
for iDate = 2:size(allDates,1)
    clear expDate datFileNum 
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];
    if ~exist([saveFolder '\IntraProbeCorr'],'dir'); [~,~] = mkdir([saveFolder '\IntraProbeCorr']); end
    if ~exist([saveFolder '\InterProbeCorr'],'dir'); [~,~] = mkdir([saveFolder '\InterProbeCorr']); end
    
    for iFile = 1:length(datFileNum)
        clear corrRefA corrRefB
        fileNum = datFileNum(iFile);
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if isempty(probeA) || isempty(probeB); continue; end
       if ~strcmp(expDate,'08_08_2022'); probeB(:,[14,22]) = []; end
       probeBList = [1:13 15:21 23:32]; 

        for iBand = 1:3
            switch iBand
                case 1
                    pA = probeA;
                    pB = probeB;
                    figTitle = 'WB';
                case 2
                    [b,a] = butter(3,(gammaBand./(fs/2)),'bandpass'); % Filtering in the gamma band
                    pA = filtfilt(b,a,probeA);
                    pB = filtfilt(b,a,probeB);
                    figTitle = 'Gamma band';
                case 3
                    [b,a] = butter(3,(alphaBand./(fs/2)),'bandpass'); % Filtering in the alpha band
                    pA = filtfilt(b,a,probeA);
                    pB = filtfilt(b,a,probeB);
                    figTitle = 'Alpha band';

            end
            clear chA chB
            chA(1) = chInCortexProbeA{iDate}(iFile);
            if chA(1)+17 > 32
                chA(2) = 32;
            else
                chA(2) = chA(1)+17;
            end

            chB(1) = chInCortexProbeB{iDate}(iFile);
            if chB(1)+17 > size(pB,2)
                chB(2) = size(pB,2);
            else
                chB(2) = chB(1)+17;
            end

            corrChA = corr(pA);
            corrChB = corr(pB);
            corrAB  = corr(pA(:,chA(1):chA(2)),pB(:,chB(1):chB(2)),'rows','complete');
            corrAll{iDate,fileNum}= corrAB;
            nanVals = isnan(corrAB); 
%             corrAB = corrAB(:,any(~isnan(corrAB))); % Remove Nan in column
%             corrAB = corrAB(any(~isnan(corrAB)),:); % Remove Nan in row

            corrRefA(:,iBand) = max(corrAB,[],2,'omitnan'); %mean(corrAB,2,'omitnan');
            corrRefB(:,iBand) =  max(corrAB,[],1,'omitnan');%mean(corrAB,1,'omitnan'); 

            % Save the intra-probe correlation heatmaps for Probe A and
            % Probe B in gamma and wideband only
            if iBand~=3
                if ~exist([saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'file')
                    figure; imagesc(imgaussfilt(corrChA,1)); title(['Probe A: Intra probe correlation -  ' figTitle]);
                    caxis([0 1]);yticks(1:32); xticks(1:32);colormap jet; colorbar;
                    f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300); close gcf;
                end

                if ~exist([saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'file')
                    figure; imagesc(imgaussfilt(corrChB,1)); title(['Probe B: Intra probe correlation -  ' figTitle]);colormap jet;
                    caxis([0 1]);yticks(1:30); xticks(1:30); xticklabels(probeBList); yticklabels([1:13 15:21 23:32]); colorbar
                    f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300); close gcf;
                end
            end

            % Save the pairwise correlation heatmaps between Probe A and B
            if ~exist([saveFolder '\InterProbeCorr\File_' num2str(fileNum) '_' figTitle '.png'],'file')
                figure; imagesc(imgaussfilt(corrAB,1)); title(['Pairwise correlation -  ' figTitle]); colormap jet;  xlabel('Probe B'); ylabel('Probe A'); colorbar;
                clim([0 1]); yticks(1:18); xticks(1:18); %xticklabels(probeBList);
                f = gcf; exportgraphics(f,[saveFolder '\InterProbeCorr\File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300); close gcf;
            end
        end

        % Save the pairwise correlation with respect to probeA or
        % probeB for different bands
%         if ~exist([saveFolder '\InterProbeCorr\RefWiseFile_' num2str(fileNum) '.png'],'file')
            figure; lineCol = {'k';'r';'b'};

            subplot(1,2,1); % Ref A
            for iB = 1:3
                plot(corrRefA(:,iB),1:size(corrRefA,1),'Color',lineCol{iB},'LineWidth',1); hold on;
            end
            set(gca, 'YDir','reverse');grid on; 
            xlabel('Correlation'); ylabel('Channels'); yticks(1:size(corrRefA,1)); title('Ref: Probe A');
            xlim([round(min(corrRefA,[],'all')-0.1,1) round(max(corrRefA,[],'all')+0.1,1)]);
            xticks(round(min(corrRefA,[],'all')-0.1:0.1:max(corrRefA,[],'all')+0.1,1));

            subplot(1,2,2); % Ref B
            for iB = 1:3
                plot(corrRefB(:,iB),1:size(corrRefB,1),'Color',lineCol{iB},'LineWidth',1); hold on;
            end
            set(gca, 'YDir','reverse'); grid on; 
            xlabel('Correlation'); ylabel('Channels'); yticks(1:size(corrRefB,1)); title('Ref: Probe B');
            xlim([round(min(corrRefB,[],'all')-0.1,1) round(max(corrRefB,[],'all')+0.1,1)]);
            xticks(round(min(corrRefB,[],'all')-0.1:0.1:max(corrRefB,[],'all')+0.1,1));

            sgtitle(strrep([expDate ' datafile ' num2str(fileNum) ],'_','\_'));  legend('WB','Gamma','Alpha','Location','Northeast');
            f = gcf; exportgraphics(f,[saveFolder '\InterProbeCorr\MaxRefWiseFile_' num2str(fileNum) '.png'],'Resolution',300); close gcf;
%         end
    end
end

%% 7. Pairwise correlation of the envelope - to get the instantaneous power
for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];
    if ~exist([saveFolder '\EnvelopeCorr'],'dir'); [~,~] = mkdir([saveFolder '\EnvelopeCorr']); end

    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};

        if isempty(probeA) || isempty(probeB); continue; end
        if ~strcmp(expDate,'08_08_2022'); probeB(:,[14,22]) = []; end
        probeBList = [1:13 15:21 23:32];
        

        for iBand = 1:3
            switch iBand
                case 1
                    pA = probeA;
                    pB = probeB;
                    figTitle = 'WB';
                case 2
                    [b,a] = butter(3,(gammaBand./(fs/2)),'bandpass'); % Filtering in the gamma band
                    pA = filtfilt(b,a,probeA);
                    pB = filtfilt(b,a,probeB);
                    figTitle = 'Gamma band';
                case 3
                    [b,a] = butter(3,(alphaBand./(fs/2)),'bandpass'); % Filtering in the alpha band
                    pA = filtfilt(b,a,probeA);
                    pB = filtfilt(b,a,probeB);
                    figTitle = 'Alpha band';

            end
            clear chA chB
            chA(1) = chInCortexProbeA{iDate}(iFile);
            if chA(1)+17 > 32
                chA(2) = 32;
            else
                chA(2) = chA(1)+17;
            end

            chB(1) = chInCortexProbeB{iDate}(iFile);
            if chB(1)+17 > size(pB,2)
                chB(2) = size(pB,2);
            else
                chB(2) = chB(1)+17;
            end
            [pA,~] = envelope(pA); [pB,~] = envelope(pB);

            % Save the pairwise correlation heatmaps between Probe A and B
            if ~exist([saveFolder '\EnvelopeCorr\Envelope_File_' num2str(fileNum) '_' figTitle '.png'],'file') || ~exist([saveFolder '\EnvelopeCorr\Envelope_File_' num2str(fileNum) '_' figTitle '.eps'],'file')
                figure; imagesc(imgaussfilt(corr(pA,pB,'rows','complete'),1)); title(['Pairwise correlation of envelope -  ' figTitle]); colormap jet;  xlabel('Probe B'); ylabel('Probe A'); colorbar;
                clim([-0.5 0.5]); yticks(1:size(pA,2)); xticks(1:size(pB,2)); %xticklabels(probeBList);
                f = gcf; exportgraphics(f,[saveFolder '\EnvelopeCorr\Envelope_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                exportgraphics(f,[saveFolder '\EnvelopeCorr\Envelope_File_' num2str(fileNum) '_' figTitle '.eps'],'ContentType','vector');
                close gcf;
            end

            corrABEnv  = corr(pA(:,chA(1):chA(2)),pB(:,chB(1):chB(2)),'rows','complete');
            corrAllEnv{iDate,fileNum}= corrABEnv;
%             corrRefA(:,iBand) = max(corrABEnv,[],2,'omitnan'); %mean(corrAB,2,'omitnan');
%             corrRefB(:,iBand) =  max(corrABEnv,[],1,'omitnan');%mean(corrAB,1,'omitnan'); 
        end
    end
end

%% Figures for testing...
for iProbe = 1:2
    switch iProbe
        case 1
            insideCortex = corrChA;
            outsideCortex = corrChAOut;
            randomCorr = corrChARandom;
            figTitle = 'Probe A';
            refCh    = chA;
        case 2
            insideCortex = corrChB;
            outsideCortex = corrChBOut;
            randomCorr = corrChBRandom;
            figTitle = 'Probe B';
            refCh    = chB;
    end
    numCh = size(insideCortex,2);

    % Plot the correlation vs channels
    meanIn = mean(squeeze(insideCortex(1,:,:)),2,'omitnan'); stdIn = std(squeeze(insideCortex(1,:,:)),[],2,'omitnan');
    figure; plot(meanIn,'LineWidth',2); hold on;
    patch([(1:numCh)' ;flipud((1:numCh)')],[(meanIn-stdIn);  flipud(meanIn+stdIn)],'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Channels'); ylabel('Correlation'); title([figTitle ' Reference channels ' num2str(refCh:refCh+2)]);
    if iProbe==2; xticks(1:30); xticklabels(probeBList); else; xticks(1:32); end

    meanOut = mean(squeeze(outsideCortex(1,:,:)),2,'omitnan'); stdOut = std(squeeze(outsideCortex(1,:,:)),[],2,'omitnan');
    figure; plot(meanOut,'LineWidth',2); hold on;
    patch([(1:numCh)' ;flipud((1:numCh)')],[(meanOut-stdOut);  flipud(meanOut+stdOut)],'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Channels'); ylabel('Correlation'); title([figTitle ' Reference channels ' num2str(chOutCortex)]);
    if iProbe==2; xticks(1:30); xticklabels(probeBList); else; xticks(1:32);end

    meanDeep = mean(squeeze(randomCorr(1,:,:)),2,'omitnan');stdDeep = std(squeeze(randomCorr(1,:,:)),[],2,'omitnan');
    figure;  plot(meanDeep,'LineWidth',2); hold on;
    patch([(1:numCh)' ;flipud((1:numCh)')],[(meanDeep-stdDeep);  flipud(meanDeep+stdDeep)],'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Channels'); ylabel('Correlation'); title([figTitle ' Reference channels ' num2str(chDeep)]);
    if iProbe==2; xticks(1:30); xticklabels(probeBList); else; xticks(1:32);end

    % Plot the slope for the channels
    slopeInside  = movmean(squeeze(gradient(insideCortex(1,:,:))),2,'omitnan');  meanSlopeIn   = mean(slopeInside,2,'omitnan');  stdSlopeIn = std(slopeInside,[],2,'omitnan');
    slopeOutside = movmean(squeeze(gradient(outsideCortex(1,:,:))),2,'omitnan'); meanSlopeOut  = mean(slopeOutside,2,'omitnan'); stdSlopeOut = std(slopeOutside,[],2,'omitnan');
    slopeDeep    = movmean(squeeze(gradient(randomCorr(1,:,:))),2,'omitnan');    meanSlopeDeep = mean(slopeDeep,2,'omitnan');    stdSlopeDeep = std(slopeDeep,[],2,'omitnan');


    figure;  plot(meanSlopeIn,'LineWidth',2); hold on;
    patch([(1:numCh)' ;flipud((1:numCh)')],[(meanSlopeIn-stdSlopeIn);  flipud(meanSlopeIn+stdSlopeIn)],'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Channels'); ylabel('Slope'); title([figTitle ' Reference channels ' num2str(refCh:refCh+2)]);
    if iProbe==2; xticks(1:30); xticklabels(probeBList);else; xticks(1:32); end

    figure;plot(meanSlopeOut,'LineWidth',2); hold on;
    patch([(1:numCh)' ;flipud((1:numCh)')],[(meanSlopeOut-stdSlopeOut);  flipud(meanSlopeOut+stdSlopeOut)],'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Channels'); ylabel('Slope'); title([figTitle ' Reference channels ' num2str(chOutCortex)]);
    if iProbe==2; xticks(1:30); xticklabels(probeBList); else; xticks(1:32);end

    figure;plot(meanSlopeDeep,'LineWidth',2); hold on;
    patch([(1:numCh)' ;flipud((1:numCh)')],[(meanSlopeDeep-stdSlopeDeep);  flipud(meanSlopeDeep+stdSlopeDeep)],'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Channels'); ylabel('Slope'); title([figTitle ' Reference channels ' num2str(chDeep)]);
    if iProbe==2; xticks(1:30); xticklabels(probeBList);else; xticks(1:32); end

end

%%
[b6,a6] = butter(3,[6 250]./(fs/2),'bandpass');
pA = probe1{fileNum,iDate};% pA = pA(1:300e3,:);
pAF = filtfilt(b6,a6, pA);
pB = probe2{fileNum,iDate};
pB(:,[14 18 22 ]) = []; probeBList = [1:13 15:21 23:32];
% pB = pB(1:300e3,:);
pBF = filtfilt(b6,a6,pB);

figure; imagesc(imgaussfilt(corr(pA(1:300e3,:)),1)); title(['Probe A: Intra probe correlation -  ' figTitle]);colormap jet;
caxis([0 1]);yticks(1:30); xticks(1:30); xticklabels([1:13 15:21 23:32]); yticklabels([1:13 15:21 23:32]); colorbar

%%
L        = size(probeA,1);
binWidth = 30e3;
w        = 300e3;
chLen    = size(pB,2);
clear corrDat
ind = 1;
for iW = 1:binWidth:L-w+1
    corrChAll  = corrcoef(pB,'rows','complete');
    corrChAllF = corrcoef(pBF,'rows','complete');

    corrDat  = corrcoef(pB(iW:iW+w-1,:),'rows','complete');
    corrDatF = corrcoef(pBF(iW:iW+w-1,:),'rows','complete');

    corrMaps(:,ind)  =  corr(reshape(squeeze(corrDat),[chLen*chLen 1]),reshape(corrChAll,[chLen*chLen 1]));
    corrMapsF(:,ind) = corr(reshape(squeeze(corrDatF),[chLen*chLen 1]),reshape(corrChAllF,[chLen*chLen 1]));

    ind = ind+1;
end
boxplot([corrMaps' corrMapsF'],'Notch','off','Labels',{'1-250 Hz'; '6-250 Hz'});
title('Correlations between 5 and 10 minute segments before and after removing 1-5 Hz')
ylabel('Correlations'); %ylim([0.95 1.01]);

%% Intra probe correlations
chA = chInCortexProbeA{iDate}(iFile);
chB = chInCortexProbeB{iDate}(iFile);

ch1 = 1:4;
ch2 = [15:17];
ch3 = 25:29;
%
% corrCh1 = corr(pB(300e3:end,:),pB(300e3:end,ch1));
% corrCh2 = corr(pB(300e3:end,:),pB(300e3:end,ch2));
% corrCh3 = corr(pB(300e3:end,:),pB(300e3:end,ch3));

corrCh1 = corr(pBF(1:300e3,:),pBF(1:300e3,ch1));
corrCh2 = corr(pBF(1:300e3,:),pBF(1:300e3,ch2));
corrCh3 = corr(pBF(1:300e3,:),pBF(1:300e3,ch3));

figure;
subplot(1,3,1); imagesc(imgaussfilt(mean(corrCh1,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:32); title('Ch:1-4'); yticklabels(probeBList);
subplot(1,3,2); imagesc(imgaussfilt(mean(corrCh2,2),1)); colormap jet; clim([0 1]); colorbar; yticks(1:32); title('Ch:15-19');yticklabels(probeBList);
subplot(1,3,3); imagesc(imgaussfilt(mean(corrCh3,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:32); title('Ch:25-29');yticklabels(probeBList);





