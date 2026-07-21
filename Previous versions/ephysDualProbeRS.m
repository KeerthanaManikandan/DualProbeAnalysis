clear;clc;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir));
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
clc;

%% 1. Initializing all the relevant parameters and flags
monkeyName   = 'CharlieSheen';
hemisphere   = 'Left';
saveDataFlag = 1;
savePlotFlag = 1;
badChProbeB  = [14 22];
fs           = 1e3; % Sampling frequency

% Load all the experiment dates, data file names for the corresponding monkey
if strcmp(monkeyName,'CharlieSheen') % Charlie Sheen
    allDates         = ['11_01_2021';'01_11_2022'; '02_28_2022'; '08_30_2022'; '10_10_2022'];
    datFileNumAll    = {1:11; 5:10; [5 6]; 1:12;1:9};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\';
    refDate      = '08_31_2021';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Charlie Sheen Combined Green 08_31_2021';
    datFileNameAll   = {'000'; 'run000'; 'datafile000'; 'datafile000'; 'datafile000'};
    chInCortexProbeA = {[14 6 6 6 7 6 6 6 6 6 1]; [13 13 13 13 13 13]; [11 11]; [10 5 8 9 5 11 11 4 4 8 8 8];   [1 5 5 7 8 9 9 10 7]};
    chInCortexProbeB = {[9 9 9 5 5 5 5 5 5 5 2]; [11 11 9 11 10 6];    [10 10]; [10 10 10 10 10 8 8 8 8 8 9 5]; [1 6 6 6 6 6 6 9 9]};
    siteConnAll      = {[0 1 1 0 1 0 0 0 1 1 0 ], [0 0 1 0 1 0],       [1 1],   [1 1 1 0 0 1 0 1 0 1 1 1],      [1 1 1 0 1 1 1 0 1 ]};    

elseif strcmp(monkeyName,'Bordeaux')
    allDates     = ['01_06_2020'; '02_10_2020']; % Bordeaux
    eegChName    = 'analog1';

elseif strcmp(monkeyName,'Whiskey') % Whiskey
    allDates         = ['08_08_2022'; '09_19_2022'; '10_17_2022']; % To edit this segment
    datFileNumAll    = {1:4; 1:11; 1:15};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\';
    refDate      = '05_09_2022';
    refImageName = 'Combined Green';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    datFileNameAll   = {'datafile_000'; 'datafile000'; 'datafile000'};
    eegChName        = 'analog1';
    chInCortexProbeA = {[10 10 10 10]; [6 6 6 6 6 9 10 10 10 10 10]; [6 1 10 10 5 5 5 5 10 10 6 7 9 8 12]};
    chInCortexProbeB = {[1 1 1 1];     [9 9 10 10 6 9 8 8 11 9 9 ];  [8 10 10 10 6 6 6 6 6 6 6 6 6 6 6 ]};
    siteConnAll      = {[1 0 0 1],     [1 0 1 1 0 1 0 1 1 0 1],      [0 1 1 1 1 1 1 1 0 1 1 1 0 1 0]};

elseif strcmp(monkeyName,'Merle') % Merle
    allDates   = '05_07_2018';
    serverPath = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\Euthanized Animal Data\371-16_Merle_SqM\Left Hemisphere\';
end


%% 2. Retrieve or store the filtered data
[probe1,probe2,eeg] = saveLFPDualProbe(monkeyName,hemisphere,allDates,datFileNameAll,datFileNumAll,serverPath);

%% 3. Get the distance and the connectivity between the sites 
[distSitesAll,connValsAll,refSiteAll, movSiteAll,greenMapRef] = getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);

%% Get the intraprobe and interprobe/pairwise correlations for the recordings
for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];
    if ~exist([saveFolder '\IntraProbeCorr'],'dir'); [~,~] = mkdir([saveFolder '\IntraProbeCorr']); end
    if ~exist([saveFolder '\InterProbeCorr'],'dir'); [~,~] = mkdir([saveFolder '\InterProbeCorr']); end

    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
%         chA = chInCortexProbeA{iDate}(iFile);
%         chB = chInCortexProbeB{iDate}(iFile);
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        if isempty(probeA) || isempty(probeB); continue; end 
        probeB(:,[14,22]) = []; 
        
        for iBand = 1:2
            switch iBand
                case 1
                    pA = probeA; 
                    pB = probeB; 
                    figTitle = 'WB';
                case 2
                    [b,a] = butter(3,([6 250]./(fs/2)),'bandpass'); % Filtering in the gamma band 
                    pA = filtfilt(b,a,probeA);
                    pB = filtfilt(b,a,probeB); 
                    figTitle = 'Gamma band';
            end
            corrChA = corr(probeA(1:300e3,:)); 
            corrChB = corr(probeB(1:300e3,:));
            corrAB  = corr(pA(1:300e3,:),pB(1:300e3,:));

            figure; imagesc(imgaussfilt(corrChA,1)); title(['Probe A: Intra probe correlation -  ' figTitle]);
            caxis([0 1]);yticks(1:32); xticks(1:32);colormap jet; colorbar
            f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300); close gcf;

            figure; imagesc(imgaussfilt(corrChB,1)); title(['Probe B: Intra probe correlation -  ' figTitle]);colormap jet; 
            caxis([0 1]);yticks(1:30); xticks(1:30); xticklabels([1:13 15:21 23:32]); yticklabels([1:13 15:21 23:32]); colorbar
            f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300); close gcf;

            figure; imagesc(imgaussfilt(corrAB,1)); title(['Pairwise correlation -  ' figTitle]); colormap jet;  xlabel('Probe B'); ylabel('Probe A'); colorbar; 
            caxis([-0.5 0.5]); yticks(1:32); xticks(1:30); xticklabels([1:13 15:21 23:32]); 
            f = gcf; exportgraphics(f,[saveFolder '\InterProbeCorr\File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300); close gcf;  
        end 
    end 
end 

%% Getting intra probe correlations
[b6,a6] = butter(3,[6 250]./(fs/2),'bandpass');
ch1 = 1:5;
ch2 = 15:20; 
ch3 = 25:30;
pA = probe1{fileNum,iDate};% pA = pA(1:300e3,:);
pAF = filtfilt(b6,a6, pA);
% pB = probe1{fileNum,iDate}; pB = pB(1:300e3,:); 

% corrCh1 = corr(pA,pA(:,ch1)); 
% corrCh2 = corr(pA,pA(:,ch2)); 
% corrCh3 = corr(pA,pA(:,ch3)); 

corrCh1 = corr(pAF,pAF(:,ch1)); 
corrCh2 = corr(pAF,pAF(:,ch2)); 
corrCh3 = corr(pAF,pAF(:,ch3));

figure; 
subplot(1,3,1); imagesc(imgaussfilt(mean(corrCh1,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:32); title('Ch:1-5');
subplot(1,3,2); imagesc(imgaussfilt(mean(corrCh2,2),1)); colormap jet; clim([0 1]); colorbar; yticks(1:32); title('Ch:15-20');
subplot(1,3,3); imagesc(imgaussfilt(mean(corrCh3,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:32); title('Ch:25-30');
%%
ch1 = [13 19:22];
ch2 = [23 24 26]; 
pA = probe1{fileNum,iDate}; pA = pA(1:300e3,:);
pAF = filtfilt(b6,a6, pA);
pB = probe2{fileNum,iDate}; pB = pB(1:300e3,:); 
pBF = filtfilt(b6,a6,pB); 

% corrCh1 = corr(pA,pA(:,ch1)); 
% corrCh2 = corr(pA,pA(:,ch2)); 
% corrCh1 = corr(pAF,pAF(:,ch1)); 
% corrCh2 = corr(pAF,pAF(:,ch2)); 


figure; 
subplot(1,2,1); imagesc(imgaussfilt(mean(corrCh1,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:32); title('Ch:13 19-22');
subplot(1,2,2); imagesc(imgaussfilt(mean(corrCh2,2),1)); colormap jet; clim([0 1]); colorbar; yticks(1:32); title('Ch:23 24 26');

%%
pB = probe2{fileNum,iDate}; pB = pB(1:300e3,:); 
pB(:,[14 22])= [];
pBF = filtfilt(b6,a6,pB); 
[b10,a10] = butter(3,[11 250]./(fs/2),'bandpass');
pBF10 = filtfilt(b10,a10,pB);

ch1 = 1:5;
ch2 = [11:16]; 
ch3 = 25:30;

% corrCh1 = corr(pB,pB(:,ch1)); 
% corrCh2 = corr(pB,pB(:,ch2)); 
% corrCh3 = corr(pB,pB(:,ch3)); 

corrCh1 = corr(pBF10,pBF10(:,ch1)); 
corrCh2 = corr(pBF10,pBF10(:,ch2)); 
corrCh3 = corr(pBF10,pBF10(:,ch3));


figure; 
subplot(1,3,1); imagesc(imgaussfilt(mean(corrCh1,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:30); title('Ch:1-5'); yticklabels([1:13 15:21 23:32]);
subplot(1,3,2); imagesc(imgaussfilt(mean(corrCh2,2),1)); colormap jet; clim([0 1]); colorbar; yticks(1:30); title('Ch:11-13 15-17'); yticklabels([1:13 15:21 23:32]);
subplot(1,3,3); imagesc(imgaussfilt(mean(corrCh3,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:30); title('Ch:25-30'); yticklabels([1:13 15:21 23:32]);



%% 4. Remove artifacts and other noise observed in the LFPs
% Plot the raw LFPs and their PSDs and spectrograms 
for iDate = 5% 1:size(allDates,1) % Get the LFP for a particular date
    clear expDate datFileNum 
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];

    if ~exist([saveFolder '\LFP_TimeSeries'],'dir'); [~,~] = mkdir([saveFolder '\LFP_TimeSeries']); end 
    if ~exist([saveFolder '\LFP_AvgSpec'],'dir'); [~,~] = mkdir([saveFolder '\LFP_AvgSpec']); end

    for iFile = 8%1:length(datFileNum) % Get the LFP for a particular recording for a particular date
        fileNum = datFileNum(iFile);
        chA = chInCortexProbeA{iDate}(iFile); 
        chB = chInCortexProbeB{iDate}(iFile); 
        for iSet = 1:2 % Get the LFP for the two probes for the recording
            clear probe
            switch iSet
                case 1
                    probe     = probe1{fileNum,iDate};
                    probeName = 'A';
                    chVal     = chA; 
                case 2
                    probe     = probe2{fileNum,iDate}; % Not removing bad electrodes 14 and 22 
                    probeName = 'B';
                    chVal     = chB; 
            end
            if isempty(probe); continue; end 

            % Get the raw LFPs
            if (~exist([saveFolder '\LFP_TimeSeries\LFP_' num2str(fileNum) probeName '.png'],'file'))
                figure;ind = 1;
                for iCh = size(probe,2):-1:1
                    plot(1:size(probe,1),probe(:,iCh)+((ind-1)*1000)); hold on;
                    ind= ind+1;
                end

                xlabel('Time(samples)'); ylabel('Channels');
                yticks(0:1000:size(probe,2)*1e3);
                yticklabels(size(probe,2):-1:1);
                xlim([0 size(probe,1)]);
                ax = gca; ax.YAxis.FontSize = 5;
                sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe ' probeName ],'_','\_'));

                f = gcf;
                exportgraphics(f,[saveFolder '\LFP_TimeSeries\LFP_' num2str(fileNum) probeName '.png'],'Resolution',300);
                close gcf;
            end

            % Get the Average power spectrum of the LFPs
            params.Fs = fs;
            params.fpass = [1 100];
            params.pad = -1;
            params.tapers = [3 5];
            params.trialave = 0; % 
            [specgram,time, freq] = mtspecgramc(probe,[5 2],params);

            if (~exist([saveFolder '\LFP_AvgSpec\LFPAvgSpec_' num2str(fileNum)  probeName '.png'],'file'))
                figure;imagesc(time',freq',10.*log10(mean(specgram,3)'));
                sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe ' probeName ],'_','\_'));
                set(gca,'YDir','normal')
                xlabel('Time (s)'); ylabel('Frequency (Hz)');
                colormap jet; shading interp; colorbar;
                caxis([ -20 70]);

                f = gcf;
                exportgraphics(f,[saveFolder '\LFP_AvgSpec\LFPAvgSpec_' num2str(fileNum)  probeName '.png'],'Resolution',300);
                close gcf;
            end

            if (~exist([saveFolder '\LFP_AvgSpec\LFPAvgSpecCortex_' num2str(fileNum)  probeName '.png'],'file'))
                figure;imagesc(time',freq',10.*log10(mean(specgram(:,:,chVal:end),3)'));
                sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe ' probeName ],'_','\_'));
                set(gca,'YDir','normal')
                xlabel('Time (s)'); ylabel('Frequency (Hz)');
                colormap jet; shading interp; colorbar;
                caxis([ -20 70]);

                f = gcf;
                exportgraphics(f,[saveFolder '\LFP_AvgSpec\LFPAvgSpecCortex_' num2str(fileNum)  probeName '.png'],'Resolution',300);
                close gcf;
            end

            if (~exist([saveFolder '\LFP_AvgSpec\LFPAvgSpecCortexOut_' num2str(fileNum)  probeName '.png'],'file'))
                figure;imagesc(time',freq',10.*log10(mean(specgram(:,:,1:chVal-1),3)'));
                sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe ' probeName ],'_','\_'));
                set(gca,'YDir','normal')
                xlabel('Time (s)'); ylabel('Frequency (Hz)');
                colormap jet; shading interp; colorbar;
                caxis([ -20 70]);

                f = gcf;
                exportgraphics(f,[saveFolder '\LFP_AvgSpec\LFPAvgSpecCortexOut_' num2str(fileNum)  probeName '.png'],'Resolution',300);
                close gcf;
            end
            
            % Get the average power spectrum 
            [psd,freqVals] = mtspectrumc(probe,params);
            if (~exist([saveFolder '\LFP_AvgSpec\AvgPSD_' num2str(fileNum)  probeName '.png'],'file'))
                figure;
                subplot(1,3,1); plot(freqVals,10.*log10(mean(psd,2)')); title('All channels'); ylim([-60 80]); xlim([0 55]); xticks(0:10:50);
                ylabel('Power (dB)'); xlabel('Frequency (Hz)');
                subplot(1,3,2); plot(freqVals,10.*log10(mean(psd(:,chVal:end),2)')); title('In cortex'); ylim([-60 80]); xlim([0 55]);xticks(0:10:50);
                subplot(1,3,3); plot(freqVals,10.*log10(mean(psd(:,1:chVal-1),2)')); title('Out of cortex'); ylim([-60 80]); xlim([0 55]);xticks(0:10:50);
                sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe ' probeName ],'_','\_'));
                f = gcf;
                exportgraphics(f,[saveFolder '\LFP_AvgSpec\AvgPSD_' num2str(fileNum)  probeName '.png'],'Resolution',300);
                close gcf;
            end
        end
    end
end
%% Comparing the spectrogram within and across the probe... 
params.Fs = fs;
params.fpass = [1 100];
params.pad = -1;
params.tapers = [3 5];
params.trialave = 0; %

probeA     = probe1{fileNum,iDate};
probeB     = probe2{fileNum,iDate};

[specgramA,timeA, freqA] = mtspecgramc(pA(1:300e3,:),[5 2],params);
[specgramB,timeB, freqB] = mtspecgramc(pB(1:300e3,:),[5 2],params);
specgramA = 10.*log10(specgramA); specgramB = 10.*log10(specgramB); % Converting to log scale...

[psdA,freqValsA] = mtspectrumc(pA(1:300e3,:),params);
[psdB,freqValsB] = mtspectrumc(pB(1:300e3,:),params);
psdA = 10.*log10(psdA); psdB = 10.*log10(psdB);

for iSet = 1:2
    clear spec p 
    switch iSet
        case 1
            spec     = specgramA; 
            p        = psdA; 
            figLabel = 'Probe A';
        case 2
            spec     = specgramB; 
            p        = psdB; 
            figTitle = 'Probe B'; 
    end
    sizeSpec = size(spec);
    corrIntraSpec(iSet,:,:) = corr(reshape(spec,[sizeSpec(1)*sizeSpec(2) sizeSpec(3)]));
    corrIntraPSD(iSet,:,:)  = corr(p); 
end

%% Plotting spectrogram and PSD of all channels

for iCh = 1:32
    figure;imagesc(timeB',freqB',specgramB(:,:,iCh)');
    sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe B - Ch ' num2str(iCh)  ],'_','\_'));
    set(gca,'YDir','normal')
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    colormap jet; shading interp; colorbar;  axis square
    caxis([ -20 70]); %ylim([1 30]);
end


%% 
chLabel = [1:13 15:21 23:32];
figure; 
for iCh = 1:32%length(chLabel)
    plot(freqValsB,psdB(:,iCh)); ylim([-60 80]); xlim([5 35]); xticks(0:0.5:50);
    xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 
    title(['Ch: ' num2str((iCh))]);% title(['Ch: ' num2str(chLabel(iCh))]);
end 
%%
chHarmonic = [11 13 20 21 23 30:32];
chNHarmonic = [12 15 16 17:19 24:29];
figure;
subplot(2,1,1); plot(freqValsB,psdB(:,chHarmonic)); ylim([-60 80]); xlim([5 35]); xticks(0:0.5:50); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title('Ch containing prominent harmonics')
subplot(2,1,2); plot(freqValsB,psdB(:,chNHarmonic)); ylim([-60 80]); xlim([5 35]); xticks(0:0.5:50); xlabel('Frequency (Hz)'); ylabel('Power (dB)');title('Ch not containing prominent harmonics')
%% 
ch1 = [11 13 20 21];
ch2 = [23 30:32]; 
% ch3 = 25:30;
% ch1 = [12 15 16 17];
% ch2 = [24:28]; 
[b10,a10] = butter(3,[11 250]./(fs/2),'bandpass');
% pN =probe2{fileNum,iDate};
pN = filtfilt(b6,a6,probe2{fileNum,iDate});

% corrCh1 = corr(psdB,psdB(:,ch1)); 
% corrCh2 = corr(psdB,psdB(:,ch2)); 
% corrCh3 = corr(psdB,psdB(:,ch3)); 

corrCh1 = corr(pN,pN(:,ch1)); 
corrCh2 = corr(pN,pN(:,ch2)); 
corrCh3 = corr(pN,pN(:,ch3));

corrCh1([14,22],:)= [];
corrCh2([14,22],:)= [];
corrCh3([14,22],:)= [];

figure; 
subplot(1,2,1); imagesc(imgaussfilt(mean(corrCh1,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:30); title('Ch:11 13 20 21'); yticklabels([1:13 15:21 23:32]);
subplot(1,2,2); imagesc(imgaussfilt(mean(corrCh2,2),1)); colormap jet; clim([0 1]); colorbar; yticks(1:30); title('Ch: 23 30:32'); yticklabels([1:13 15:21 23:32]);
% subplot(1,3,3); imagesc(imgaussfilt(mean(corrCh3,2),1)); colormap jet; clim([0 1]); colorbar;yticks(1:30); title('Ch:25-30'); yticklabels([1:13 15:21 23:32]);

%% 
figure; 

for iCh = 1:32
    clf
    plot(freqValsA,psdA(:,iCh)); ylim([-60 80]); xlim([ 1.2 3.2]); xticks(0:0.5:50);
    title(['Ch: ' num2str(iCh)]);
%     if iCh == 29
    xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 
%     end
end 

%% 
% bw = (3/(fs/2))/35;
% [bL,aL] = iircomb(fs/4,bw,'notch');
[bL,aL] = butter(3,[3 4]./(fs/2),'stop');
probeT = filtfilt(bL,aL,probe);
[psd,freqVals] = mtspectrumc(probeT,params);
[specgram,time, freq] = mtspecgramc(probeT,[5 2],params);


figure;
subplot(1,3,1); plot(freqVals,10.*log10(mean(psd,2)')); title('All channels'); ylim([-60 80]); xlim([0 55]); xticks(0:10:50);
subplot(1,3,2); plot(freqVals,10.*log10(mean(psd(:,chVal:end),2)')); title('In cortex'); ylim([-60 80]); xlim([0 55]);xticks(0:10:50);
subplot(1,3,3); plot(freqVals,10.*log10(mean(psd(:,1:chVal-1),2)')); title('Out of cortex'); ylim([-60 80]); xlim([0 55]);xticks(0:10:50);

%% Remove the time segment 125:150 s; 375:450 s

probeN = probe([1:125e3 150e3:375e3 450e3:end],:);
[specgramN,timeN, freqN] = mtspecgramc(probe,[5 2],params); 
figure;imagesc(timeN',freqN',10.*log10(specgramN(:,:,18)'));

set(gca,'YDir','normal')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
colormap jet; shading interp; colorbar
caxis([ -20 70]);

[psdN,freqValsN] = mtspectrumc(probeN,params);
figure; plot(freqValsN,10.*log10(psdN(:,18)));
xlim([ 0 55]); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 

%% Getting the average spectrogram 
[specgramN,timeN, freqN] = mtspecgramc(probe,[5 2],params); 
figure;imagesc(timeN',freqN',10.*log10(mean(specgramN,3)'));
set(gca,'YDir','normal')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
colormap jet; shading interp; colorbar
caxis([ -20 70]);
title('Average of all 32 channels'); 

%% Removing 400-500 second segment 
probeT = probe([1:510e3 520e3:end],:);
probeT(:,badChProbeB) = [];
% bw = (10/(fs/2))/35;
% % [bL,aL] = butter(3,[9.5 10.5]./(fs/2),'stop');
% [bL,aL] = iircomb(fs/10,bw,'notch'); probeT = filtfilt(bL,aL,probeT); 

[specgramT,timeT, freqT] = mtspecgramc(probeT,[5 2],params); 
figure;imagesc(timeT',freqT',10.*log10(mean(specgramT,3)'));
set(gca,'YDir','normal')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
colormap jet; shading interp; colorbar
caxis([ -20 70]);
title('Average of all 32 channels'); 

%%
probeT = probe([1:360e3 500e3:end],:);
probeT(:,badChProbeB) = []; bw = (10/(fs/2))/35;
% % [bL,aL] = butter(3,[9.5 10.5]./(fs/2),'stop');
% [bL,aL] = iircomb(fs/10,bw,'notch'); 
% probeT = filtfilt(bL,aL,probeT); 
figure;ind = 1;
for iCh = size(probeT,2):-1:1
    plot(1:size(probeT,1),probeT(:,iCh)+((ind-1)*1000)); hold on;
    ind= ind+1;
end
xlabel('Time(samples)'); ylabel('Channels');
yticks(0:1000:size(probeT,2)*1e3);
yticklabels(size(probeT,2):-1:1);
xlim([0 size(probeT,1)]);

sgtitle(strrep([expDate ' datafile:' num2str(fileNum) ' Probe ' probeName ],'_','\_'));


%% Get the intra-probe correlations
chOutCortex = 1:2;
chRandom     = [ 10 12];
chA          = chInCortexProbeB{5}(2)+10; % 1mm from the transition point listed in experiment notes
gammaBand    = [30 90];
[b,a]        = butter(3,gammaBand./(fs/2),'bandpass'); 

for iF = 1:2
    switch iF
        case 1
            filterFlag = 0; 
            bandName   = 'WB';
            probeTemp  = probeT;
%             probeTemp(:,badChProbeB) =[];
        case 2
            filterFlag = 1;
            bandName   = 'Gamma';
            probeTemp  = probeT;
%             probeTemp(:,badChProbeB) =[];
            probeTemp  = filtfilt(b,a,probeTemp);
    end
    cortexIn     = probeTemp(:,chA:chA+1);               % Correlating with channels that we know is inside the cortex
    cortexOut    = probeTemp(:,chOutCortex);             % Correlating with channels that we know is out of cortex
    cortexRandom = probeTemp(:,chRandom);                % Correlating with random channels

    corrCh       = corr(probeTemp,cortexIn);
    corrChOut    = corr(probeTemp,cortexOut);
    corrChRandom = corr(probeTemp,cortexRandom); % Self = 1; across 2;

    figure;
    subplot(1,3,1); imagesc(movmean(mean(corrCh,2,'omitnan'),5));  caxis([ 0 1]); colormap jet;
    ylabel('Mean correlation'); title('Inside cortex ');     yticks(1:size(probeT,2)); colorbar;

    subplot(1,3,2); imagesc(movmean(mean(corrChOut,2,'omitnan'),5)); caxis([ 0 1]); colormap jet;
    ylabel('Mean correlation');title( 'Outside cortex '); yticks(1:size(probeT,2)); colorbar;

    subplot(1,3,3); imagesc(movmean(mean(corrChRandom,2,'omitnan'),5));    caxis([0 1]); colormap jet;
    ylabel('Mean correlation'); title('Random channels'); yticks(1:size(probeT,2)); colorbar;

    sgtitle(strrep([expDate ' Datafile ' num2str(iFile)  ' - ' bandName  ],'_','\_'));
end


%%  Bad channel/recording
for iCh = 1:size(probe,2)
    corrSpec(:,iCh) = mean(corr(specgram(:,:,1),specgram(:,:,iCh)),'all'); 
end
meanVal = mean(corrSpec(3:end),'omitnan');
stdVal  = std(corrSpec(3:end),'omitnan');
figure; plot(corrSpec,'*-','MarkerSize',8); hold on; 
yline(meanVal+stdVal,'--','LineWidth',1); 
yline(meanVal-stdVal,'--','LineWidth',1);
xticks(1:32); ylabel('Correlations'); xlabel('Channels');

%% Noise affecting all channels

meanProbe = mean(probe,2,'omitnan'); 
stdProbe  = std(probe,[],2,'omitnan');

figure; plot(probe(:,18)); hold on; yline(meanProbe+stdProbe,'--','LineWidth',1); 
yline(meanProbe-stdProbe,'--','LineWidth',1);

%% Noise specific to a particular channel
t = 120e3; % 120 seconds
binWidth = 60e3; % 60 seconds
len = size(probe,1);
for iCh = 18%1:size(probe,2)
    chSpec = probe(:,iCh);%(:,:,iCh);
    ind = 1;
    binVals = 1:binWidth:(len-t+1);
    for iW = 1:length(binVals)
        if iW>1
            prevVal = binVals(iW-1);
            val     = binVals(iW);
                corrSpec(iW,:) = mean(corr(chSpec(prevVal:prevVal+t-1,:),chSpec(val:val+t-1,:)),'all');
        else
            corrSpec(iW,:) = 1;
        end     
    end     
end 
figure; plot(binVals,corrSpec,'*-','MarkerSize',5);

%% 5. Find which channels are inside vs outside cortex
% Obtain the intra (within) and inter (across)  probe correlation for a set
% of channels that are known to be within cortex for wideband only
% Plot the correlation and find the cliff - either using a threshold or by
% fitting the curve and finding the transition point
% Get 18 contacts after finding the transition point - thickness of cortex
% is assumed  to be 1800 microns or 1.8mm
clear b a
filterFlag = 0;
gammaBand = [30 90]; [b,a] = butter(3,gammaBand./(fs/2),'bandpass');

for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    %     if strcmp(expDate,'08_08_2022'); continue; end

    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        clear corrChA corrChAOut corrChARandom corrChB corrChBOut corrChBRandom
        if strcmp(expDate,'09_19_2022') && fileNum == 4; continue; end    
        probeA = probe1{fileNum,iDate};
        probeB = probe2{fileNum,iDate};
        if isempty(probeA) && isempty(probeB); continue; end 

         if ~strcmp(expDate,'08_08_2022'); probeB(:,badChProbeB)= []; end

        if filterFlag
            probeA   = filtfilt(b,a,probeA);
            probeB   = filtfilt(b,a,probeB);
            bandName = 'Gamma';
        else
            bandName = 'WB';
        end
        chOutCortex = 1:3;
        chRandom    = [ 5 8 29];

        chA       = chInCortexProbeA{iDate}(iFile)+10; % 1mm from the transition point listed in experiment notes        
        cortexChA = probeA(:,chA:chA+2);               % Correlating with channels that we know is inside the cortex     
        chOutA    = probeA(:,chOutCortex);             % Correlating with channels that we know is out of cortex
        chRandomA = probeA(:,chRandom);                % Correlating with random channels

        if ~strcmp(expDate,'08_08_2022')
            chB       = chInCortexProbeB{iDate}(iFile)+10;
            cortexChB = probeB(:,chB:chB+2);
            chOutB    = probeB(:,chOutCortex);
            chRandomB = probeB(:,chRandom);
        else
            cortexChB = probeB;
            chOutB    = zeros(size(probeB));
            chRandomB = probeB;
        end

        clear fxA fxB 
        for iSet = 1:2 % Within and across probes correlation
            clear cortexIn cortexOut cortexRandom
            switch iSet
                case 1
                    cortexIn = cortexChA;
                    cortexOut = chOutA;
                    cortexRandom = chRandomA;

                case 2
                    cortexIn = cortexChB;
                    cortexOut = chOutB;
                    cortexRandom = chRandomB;
            end

            if iSet == 1; ind = 2; else ind = 1; end 

            if ~(size(cortexIn,2)== 1)
                corrChA(iSet,:,:) = corr(probeA,cortexIn);
                corrChAOut(iSet,:,:) = corr(probeA,cortexOut);
                corrChARandom(iSet,:,:) = corr(probeA,cortexRandom); % Self = 1; across 2;
            end

            if ~(size(probeB,2) == 1)
                corrChB(ind,:,:) = corr(probeB,cortexIn);
                corrChBOut(ind,:,:) = corr(probeB,cortexOut);
                corrChBRandom(ind,:,:) = corr(probeB,cortexRandom);
            else
                corrChB(ind,:,:) = corr(probeB,cortexIn)';
                corrChBOut(ind,:,:) = corr(probeB,cortexOut)';
                corrChBRandom(ind,:,:) = corr(probeB,cortexRandom)';
            end
        end

        % Taking the within correlations and computing the slope, only
        % considering channels that are known to be inside and out of
        % cortex
        fxA = movmean(mean(squeeze(gradient(corrChA(1,:,:))),2),2); [~,maxA] = max(fxA(1:15));
        fxAOut= movmean(mean(squeeze(gradient(corrChAOut(1,:,:))),2),2); [~,maxAOut] = min(fxAOut(1:15));

        if ~(size(probeB,2) == 1)
            fxB = movmean(mean(squeeze(gradient(corrChB(1,:,:))),2),2); [~,maxB] = max(fxB(1:15));
            fxBOut = movmean(mean(squeeze(gradient(corrChBOut(1,:,:))),2),2);[~,maxBOut] = min(fxBOut(1:15));
        else
            maxB = 1;
            maxBOut = 1;  
        end
        if ~filterFlag
            % Obtaining the transition point
            % If the difference between maxA and maxAOut are greater than or equal to 3, minimum of the two is considered to be the transition
            estChInCortexA{iDate}(fileNum,1) = round(mean([maxA maxAOut]));

            if estChInCortexA{iDate}(fileNum,1)+18 >size(probeA,2)
                estChInCortexA{iDate}(fileNum,2) = size(probeA,2);
            else
                estChInCortexA{iDate}(fileNum,2) = estChInCortexA{iDate}(fileNum,1)+18;
            end

            % If the difference between maxB and maxBOut are greater than or equal to 3, minimum of the two is considered to be the transition
            if ~strcmp(expDate,'08_08_2022')
                estChInCortexB{iDate}(fileNum,1) = round(mean([maxB maxBOut]));

                if estChInCortexB{iDate}(fileNum,1)+18 >size(probeB,2)
                    estChInCortexB{iDate}(fileNum,2) = size(probeB,2);
                else
                    estChInCortexB{iDate}(fileNum,2) = estChInCortexB{iDate}(fileNum,1)+18;
                end
            else
                estChInCortexB{iDate}(fileNum,1) = 1; estChInCortexB{iDate}(fileNum,2) = 1;
            end
        end
        for iProbe = 1:2
            if size(probeB,2)== 1 && iProbe == 2; continue; end
            for iSet = 1:2
                clear insideCortex outsideCortex randomCorr figTitle
                figure;
                switch iProbe
                    case 1
                        insideCortex = corrChA;
                        outsideCortex = corrChAOut;
                        randomCorr = corrChARandom;
                        figTitle = 'Probe A';
                    case 2
                        insideCortex = corrChB;
                        outsideCortex = corrChBOut;
                        randomCorr = corrChBRandom;
                        figTitle = 'Probe B';
                end

                if iSet == 1; typeTitle = 'Within'; else; typeTitle = 'Across'; end

                if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'],'dir')
                    [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
                end

                if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\Transition_' num2str(fileNum) '_' figTitle '_' typeTitle ' ' bandName '.png'],'file')
                    try mean(squeeze(insideCortex(iSet,:,:)));
                        subplot(1,3,1); imagesc(movmean(mean(squeeze((insideCortex(iSet,:,:))),2,'omitnan'),5));  caxis([ 0 1]); colormap jet;
                        ylabel('Mean correlation'); title('Inside cortex ');     yticks(1:size(insideCortex,2)); colorbar;

                        subplot(1,3,2); imagesc(movmean(mean(squeeze((outsideCortex(iSet,:,:))),2,'omitnan'),5)); caxis([ 0 1]); colormap jet;
                        ylabel('Mean correlation');title( 'Outside cortex '); yticks(1:size(insideCortex,2)); colorbar;

                        subplot(1,3,3); imagesc(movmean(mean(squeeze((randomCorr(iSet,:,:))),2,'omitnan'),5));    caxis([0 1]); colormap jet;
                        ylabel('Mean correlation'); title('Random channels'); yticks(1:size(insideCortex,2)); colorbar;

                        sgtitle(strrep([expDate ' Datafile ' num2str(iFile) ' ' figTitle ' - ' bandName ' ' typeTitle ],'_','\_'));

                        f = gcf;
                        exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\Transition_' num2str(fileNum) '_' figTitle '_' typeTitle ' ' bandName '.png'],'Resolution',300);
                        close gcf;
                    catch
                        close gcf; 
                        continue;
                    end
                else
                    close gcf;
                end
            end
        end
    end
end

%% Save the correlogram heatmap of the dual probes
filterFlag = 0; 
for iDate = 1:size(allDates,1)
    datFileNum = datFileNumAll{iDate,1};
    expDate = allDates(iDate,:);
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        for iSet = 1:2
            clear probe pow freq figTitle chInCortex
            switch iSet
                case 1
                    probe = probe1{fileNum,iDate};
                    figTitle = 'Probe A';
                    chInCortex = chInCortexProbeA{iDate};
                case 2
                    probe = probe2{fileNum,iDate};
                    if size(probe,2)>1 % To separate the single channel probe
                        probe(:,badChProbeB)= [];
                    end
                    figTitle = 'Probe B';
                    chInCortex = chInCortexProbeB{iDate};
            end
            if filterFlag % Filtering gamma band
                probe = filtfilt(b,a,probe);
                bandName = 'Gamma';
            else
                bandName = 'WB'; % Wideband
            end
            if isempty(probe); continue; end 

            corrChannels = corrcoef(probe,'rows','complete');
            figure; imagesc(imgaussfilt(corrChannels,1)); colormap jet; axis image; colorbar; hold on;
            caxis([0 1]); xticks(1:length(channels)); yticks(1:length(channels));
            sgtitle(strrep([ 'Datafile: ' num2str(fileNum) ' ' figTitle],'_','\_'));

            f = gcf;
            if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
                [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
            end
            exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\' bandName '_heatmap_' num2str(fileNum) '_' num2str(iSet) '.png'],'Resolution',300);
            close gcf;

        end
    end
end

%% Correlation vs Distance vs Connectivity 
for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:); 
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        if (strcmp(expDate,'09_19_2022') && fileNum == 4) || strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7) % Datafile 4 from 09/19/2022
            corrProbes{fileNum,iDate}     = NaN;
            meanCorrProbes(fileNum,iDate) = NaN;
            continue;
        end
        if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate}); continue; end

        clear xA yA
        xA = probe1{fileNum,iDate}(:,estChInCortexA{iDate}(fileNum,1):estChInCortexA{iDate}(fileNum,2));
        if chInCortexProbeB{iDate}(iFile)~= 1
            yA = probe2{fileNum,iDate};
            yA(:,badChProbeB)= [];
            yA = yA(:,estChInCortexB{iDate}(fileNum,1):estChInCortexB{iDate}(fileNum,2));
        else
            yA = probe2{iFile,iDate};
        end

        probe1Corr(fileNum,iDate) = mean(corr(xA),[1 2],'omitnan'); 
        probe2Corr(fileNum,iDate) = mean(corr(yA),[1 2],'omitnan');

        % Append zeros if number of channels inside the cortex is not equal
        if size(xA,2)>size(yA,2)
            yA = [yA NaN(size(yA,1), size(xA,2)-size(yA,2))];
        else
            xA = [xA NaN(size(xA,1), size(yA,2)-size(xA,2))];
        end
        corrProbes{fileNum,iDate} = corr(xA,yA);
        meanCorrProbes(fileNum,iDate) = mean(corrProbes{fileNum,iDate},[1 2],'omitnan'); 
    end
end
probe1Corr(probe1Corr==0) = []; probe2Corr(probe2Corr==0) = [];
probe1Corr = probe1Corr'; probe2Corr = probe2Corr'; 
% 2D plots of Correlation vs Distance and Correlation vs Connectivity - To
% add Area as another dimension
for iVal = 1:2
    clear valToPlot figTitle xVal
    switch iVal
        case 1
            valToPlot = cell2mat(distSitesAll);
            figTitle = 'Correlation vs Distance';
            xVal = 'Distance (mm)';
        case 2
            valToPlot = cell2mat(connValsAll);
            figTitle = 'Correlation vs Connectivity';
            xVal = 'Connectivity';
    end
    siteConnAllTemp = logical(cell2mat(siteConnAll))'; 
    if strcmp(expDate,'11_01_2021'); siteConnAllTemp(2)= []; end % Removing the 2nd datafile since that was not read by neuroshare
    meanCorrProbesTemp = meanCorrProbes; meanCorrProbesTemp(meanCorrProbesTemp==0)= []; meanCorrProbesTemp = meanCorrProbesTemp';
    nanIdx = find(isnan(meanCorrProbesTemp)); meanCorrProbesTemp(nanIdx) =[]; siteConnAllTemp(nanIdx) = []; valToPlot(nanIdx) = []; 
    figure;
    for iPlot = 1:2
        clear allVals corrVals colorVals
        switch iPlot
            case 1
                if iVal == 1
                    allVals =  [zeros(size(probe1Corr)); valToPlot(siteConnAllTemp)];
                    corrVals = [probe1Corr; meanCorrProbesTemp(siteConnAllTemp)];
                    modelType = 'exp2';
                else
                    allVals =  valToPlot(siteConnAllTemp);
                    corrVals =  meanCorrProbesTemp(siteConnAllTemp);
                    modelType = 'poly1'; axesLim = [-0.3 0.6];
                end
                colorVals =  [ 1 0 0 ];
            case 2
                if iVal == 1
                    allVals =  [zeros(size(probe1Corr)); valToPlot(~siteConnAllTemp)];
                    corrVals = [probe1Corr; meanCorrProbesTemp(~siteConnAllTemp)];
                    modelType = 'exp2';
                else
                    allVals =  valToPlot(~siteConnAllTemp);
                    corrVals =  meanCorrProbesTemp(~siteConnAllTemp);
                    modelType = 'poly1';
                end
                colorVals = [0 0 1];
        end
            plot(allVals,corrVals,'*','Color',colorVals,'MarkerSize',10);grid on; hold on;
            xFit = linspace(min(allVals),max(allVals),1000);
            [f,gof]  = fit(allVals,corrVals,modelType);
            c = coeffvalues(f); r2 = gof.rsquare;
            if iVal == 1
                fitLine = c(1)*exp(c(2)*xFit) + c(3)*exp(c(4)*xFit);
            else
                fitLine = c(1)*xFit + c(2); 
            end   
            plot(xFit,fitLine,'Color',colorVals,'LineStyle','--','LineWidth',1); 
    end
    xlabel(xVal); ylabel('Correlation value'); % legend({'Connected';'Not-connected'},'Location','northeast');
    if iVal == 2; xlim(axesLim); ylim(axesLim);else; ylim([-0.3 1]);  end 
end

%% 3D scatter of correlation, connectivity and distance
figure; 
for iDate = 1:size(allDates,1)
    clear datFileNum
    datFileNum = datFileNumAll{iDate,1};
    for iPlot = 1:2
        clear allDist corrVals colorVals connVals
        switch iPlot
            case 1
                allDist =  distSitesAll{iDate}(logical(siteConnAll{iDate}));
                corrVals = meanCorrProbes((logical(siteConnAll{iDate})),iDate);
                connVals = connValsAll{iDate}(logical(siteConnAll{iDate}));
                colorVals =  [ 1 0 0 ];
            case 2
                allDist =  distSitesAll{iDate}(~logical(siteConnAll{iDate}));
                corrVals = meanCorrProbes((~logical(siteConnAll{iDate})),iDate);
                connVals = connValsAll{iDate}(~logical(siteConnAll{iDate}));
                colorVals = [0 0 1];
        end
        plot3(allDist,corrVals,connVals,'*','Color',colorVals,'MarkerSize',10); hold on;
%         f = fit([allDist,corrVals],connVals,'poly22');
    end
end
xlabel('Distance (mm)'); ylabel('Correlation'); zlabel('Connectivity values')
legend({'Connected';'Not-connected'},'Location','northwest'); grid on;

%% Correlation values for different frequency bands
clear b a corrProbesBandWise meanCorrProbesBandWise connBands notConnBands
alphaBand = [8 12];
betaBand = [15 30];
gammaBand = [30 90];
connBands = []; notConnBands = [];

% Get the mean correlation values for different frequency bands
for iDate = 1:size(allDates,1)
    clear datFileNum expDate
    expDate = allDates(iDate,:); 
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        if strcmp(expDate,'09_19_2022') && fileNum == 4; continue; end
        if strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7); continue; end
        clear xTA yTA
        xTA = probe1{fileNum,iDate}(:,estChInCortexA{iDate}(fileNum,1):estChInCortexA{iDate}(fileNum,2));

        if chInCortexProbeB{iDate}(fileNum)~= 1
            yTA = probe2{fileNum,iDate};
            yTA(:,badChProbeB)= [];
            yTA = yTA(:,estChInCortexB{iDate}(fileNum,1):estChInCortexB{iDate}(fileNum,2));
        else
            yTA = probe2{fileNum,iDate};
        end

        % Get the filter parameters
        for iBand = 1:4
            clear b a
            switch iBand
                case 1
                    [b,a] = butter(3,alphaBand./(fs/2),'bandpass');
                case 2
                    [b,a] = butter(3,betaBand./(fs/2),'bandpass');
                case 3
                    [b,a] = butter(3,gammaBand./(fs/2),'bandpass');
            end

            % Filter them based on different frequency bands
            clear xA yA
            if iBand~=4
                xA = filtfilt(b,a,xTA);
                yA = filtfilt(b,a,yTA);
            else
                xA = xTA; yA = yTA;
            end

            % Append zeros if number of channels inside the cortex is not equal
            if size(xA,2)>size(yA,2)
                yA = [yA NaN(size(yA,1), size(xA,2)-size(yA,2))];
            else
                xA = [xA NaN(size(xA,1), size(yA,2)-size(xA,2))];
            end

            % Get the correlation values
            corrProbesBandWise{iDate,fileNum}(:,:,iBand) = corr(xA,yA);
            meanCorrProbesBandWise(iDate,fileNum,iBand) = squeeze(mean(corrProbesBandWise{iDate,fileNum}(:,:,iBand),[1 2],'omitnan'));
        end

        if strcmp(expDate,'09_19_2022') && fileNum == 4
            meanCorrProbesBandWise(iDate,fileNum,:) = NaN;
            corrProbesBandWise{iDate,fileNum} = NaN;
        end
    end

    connBands    = [connBands; squeeze(meanCorrProbesBandWise(iDate,logical(siteConnAll{iDate}),:))];
    notConnBands = [notConnBands; squeeze(meanCorrProbesBandWise(iDate,~logical(siteConnAll{iDate}),:))];
end


% Summary plot
figure; boxplot([connBands [notConnBands; NaN(size(connBands,1)-size(notConnBands,1),size(notConnBands,2))]],'Notch','off',...
    'Labels',{'Alpha band -C', 'Beta band -C','Gamma band -C','Wideband -C','Alpha band -NC', 'Beta band -NC','Gamma band -NC','Wideband -NC'});
ylabel('Correlation values');

% Sort out the mean correlation plots bandwise
alphaAll    = [connBands(:,1) [notConnBands(:,1); NaN(size(connBands,1)-size(notConnBands,1),1)]];
betaAll     =  [connBands(:,2) [notConnBands(:,2); NaN(size(connBands,1)-size(notConnBands,1),1)]];
gammaAll    =  [connBands(:,3) [notConnBands(:,3); NaN(size(connBands,1)-size(notConnBands,1),1)]];
widebandAll =  [connBands(:,4) [notConnBands(:,4); NaN(size(connBands,1)-size(notConnBands,1),1)]];

figure; boxplot([alphaAll betaAll gammaAll widebandAll ],'Notch','off',...
    'Labels',{'Alpha band -C', 'Alpha band -NC','Beta band -C', 'Beta band -NC','Gamma band -C','Gamma band -NC','Wideband -C','Wideband -NC'});
ylabel('Correlation values');

%% Pairwise correlations between probes
% Get the mean/median pairwise correlation values for different frequencies
% for different channels on the reference probe

for iDate = 1:size(allDates,1)
    clear datFileNum expDate
    expDate = allDates(iDate,:); 
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        if strcmp(expDate,'09_19_2022') && fileNum == 4; continue; end
        if strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7); continue; end
        clear xTA yTA yTB xTB
        xTA = probe1{fileNum,iDate};  %probe1{iFile,iDate}(:,estChInCortexA{iDate}(iFile,1):estChInCortexA{iDate}(iFile,2));
        xTB = probe1{fileNum,iDate}(:,estChInCortexA{iDate}(fileNum,1)+10 : estChInCortexA{iDate}(fileNum,1)+12 );

        if chInCortexProbeB{iDate}(iFile)~= 1
            yTB = probe2{fileNum,iDate};
            yTB(:,badChProbeB) = [];
            yTA = yTB(:,estChInCortexB{iDate}(fileNum,1)+10 : estChInCortexB{iDate}(fileNum,1)+12 );    
        else
            yTA = probe2{fileNum,iDate};
            yTB = probe2{fileNum,iDate}; 
        end

        % Get the filter parameters
        for iBand = 1:4
            clear b a
            switch iBand
                case 1
                    [b,a] = butter(3,alphaBand./(fs/2),'bandpass');
                case 2
                    [b,a] = butter(3,betaBand./(fs/2),'bandpass');
                case 3
                    [b,a] = butter(3,gammaBand./(fs/2),'bandpass');
            end

            % Filter them based on different frequency bands
            clear xA yA xB yB
            if iBand~=4
                xA = filtfilt(b,a,xTA);
                yA = filtfilt(b,a,yTA);
                xB = filtfilt(b,a,xTB);
                yB = filtfilt(b,a,yTB); 
            else
                xA = xTA; yA = yTA;
                xB = xTB; yB = yTB; 
            end

            % Append NaNs if number of channels inside the cortex is not equal
            if size(xA,2)>size(yA,2)
                yA = [yA NaN(size(yA,1), size(xA,2)-size(yA,2))];
            else
                xA = [xA NaN(size(xA,1), size(yA,2)-size(xA,2))];
            end
            if size(xB,2)>size(yB,2)
                yB = [yB NaN(size(yB,1), size(xB,2)-size(yB,2))];
            else
                xB = [xB NaN(size(xB,1), size(yB,2)-size(xB,2))];
            end
            % Get the correlation values
            pairWiseCorrAB{iDate,fileNum}(:,:,iBand) = corr(xA,yA);
            refProbeACorr{iDate,fileNum}(:,iBand) = mean(pairWiseCorrAB{iDate,fileNum}(:,:,iBand),2,'omitnan');

            pairWiseCorrBA{iDate,fileNum}(:,:,iBand) = corr(yB,xB); 
            refProbeBCorr{iDate,fileNum}(:,iBand) = mean(pairWiseCorrBA{iDate,fileNum}(:,:,iBand),2,'omitnan'); 
        end

        if strcmp(expDate,'09_19_2022') && fileNum == 4
            refProbeACorr{iDate,fileNum}    = NaN;
            pairWiseCorrAB{iDate,fileNum} = NaN;
            pairWiseCorrBA{iDate,fileNum}(:,:,iBand) = NaN;
            refProbeBCorr{iDate,fileNum}(:,iBand) = NaN;

        end
    end
%     connBands    = [connBands; squeeze(meanCorrProbesBandWise(iDate,logical(siteConnAll{iDate}),:))];
%     notConnBands = [notConnBands; squeeze(meanCorrProbesBandWise(iDate,~logical(siteConnAll{iDate}),:))];
end
%% Plotting 
for iDate = 2%1:size(allDates,1)
    clear datFileNum expDate
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iFile = 11% 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        if strcmp(expDate,'09_19_2022') && fileNum == 4; continue; end
        if strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 || fileNum == 7); continue; end
        for iRef = 1:2
            clear refProbeCorr chCortex figTitle
            switch iRef
                case 1
                    refProbeCorr = refProbeACorr;
                    chCortex     = estChInCortexA; 
                    figTitle     = 'ProbeA Ref'; 
                case 2
                    refProbeCorr = refProbeBCorr;
                    chCortex     = estChInCortexB; 
                    figTitle     = 'ProbeB Ref';
            end 

            figure;
            for iBand = 1:4
                switch iBand
                    case 1
                        bandName = 'Alpha';
                    case 2
                        bandName = 'Beta';
                    case 3
                        bandName = 'Gamma';
                    case 4
                        bandName = 'Wideband';
                end

                subplot(2,2,iBand);
                plot(refProbeCorr{iDate,fileNum}(:,iBand),1:size(refProbeCorr{iDate,fileNum},1));
                set(gca,'YDir','reverse');
                yline(chCortex{iDate}(fileNum,1),'LineWidth',1.5);
                yline(chCortex{iDate}(fileNum,2),'LineWidth',1.5);
                xlim([ -0.5 1]);
                xlabel('Pairwise correlation'); ylabel('Ref probe channels');
                yticks(1:4:size(refProbeCorr{iDate,fileNum},1));
                title(bandName);
            end
            sgtitle(strrep([expDate ' ' 'datafile: ' num2str(fileNum) ' '  figTitle],'_','\_'));
%             f = gcf;
%             if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'],'dir');[~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']); end
%             exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\PairwiseCorr_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
%             close gcf;
        end
    end
end

%% Minimum duration of the signal that is required for recording
% Testing with gamma band and wideband ranges
gammaBand  = [30 90];
[b,a]      = butter(3,[30 90]./(fs/2),'bandpass');
binWidth   = 60*fs;  % 1 minute bin width
filterFlag = 0;
saveFlag   = 0;

for iDate = 2%1:size(allDates,1) % Get all experiment dates
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
if strcmp(expDate,'08_08_2022'); continue; end 
    for iFile = 11%1:length(datFileNum) % Get all datafiles for each experiment date
        fileNum = datFileNum(iFile);
        if strcmp(expDate,'09_19_2022') && fileNum == 4; continue; end
        for iSet = 1:2
            clear probe pow freq figTitle chInCortex corrChAll tVals L chLen
            switch iSet
                case 1
                    probe = probe1{fileNum,iDate}(:,estChInCortexA{iDate}(fileNum,1):estChInCortexA{iDate}(fileNum,2));
                    figTitle = 'Probe A';
                    chInCortex = chInCortexProbeA{iDate};
                    figure;
                case 2
                    probe = probe2{fileNum,iDate};
                    if size(probe,2)>1 % To separate the single channel probe
                        chVals = estChInCortexB{iDate}(fileNum,1):estChInCortexB{iDate}(fileNum,2);
                        chVals(ismember(chVals,badChProbeB)) = [];                   
                        probe = probe(:,chVals);
                       
                    end
                    figTitle = 'Probe B';
                    chInCortex = chInCortexProbeB{iDate};
            end
            if filterFlag
                probe = filtfilt(b,a,probe);
                bandName = 'Gamma';
            else
                bandName = 'WB'; % Wideband
            end

            % Get the 15 minute correlation heatmap
            corrChAll = corrcoef(probe,'rows','complete');
            tVals     = (1:15).*60*fs; % Converting time to samples
            L         = size(probe,1); % Length of the signal collected
            chLen     = size(probe,2); % # of channels present in a recording
            if sum(tVals>L)>=1; tVals(tVals>L)=[]; end

            % Get the correlation heatmap for the recordings in a sliding
            % window with sliding window bin of 1 minute and variable bin
            % width
            clear corrMaps winAvgCorr
            for iT = 1:length(tVals)
                w = tVals(iT);
                clear corrDat
                ind = 1;
                for iW = 1:binWidth:L-w+1
                    if ~(size(probe,2)== 1)
                        corrDat(ind,:,:) = corrcoef(probe(iW:iW+w-1,:),'rows','complete');
                        corrMaps{iT,1}(:,ind) = corr(reshape(squeeze(corrDat(ind,:,:)),[chLen*chLen 1]),reshape(corrChAll,[chLen*chLen 1]));

                    else % Single channel condition
                        corrDat(ind,:,:) = corrcoef(probe(iW:iW+w-1,:),'rows','complete');
                        corrMaps{iT,1}(:,ind) = 1;
                    end
                    ind = ind+1;
                end
                winAvgCorr{iT,1} = mean(corrDat,[2,3],'omitnan');
            end

            corrMapAvg = cellfun(@mean,corrMaps);
            corrMapStd = cellfun(@std,corrMaps);
            yLabel     = 'Correlation';

            subplot(1,2,iSet) ;
            for iL = 1:length(corrMaps)
                plot(iL.*ones(size(corrMaps{iL,1})),corrMaps{iL,1},'b.','MarkerSize',10); hold on;
                ylim([ min(corrMapAvg)-0.2 max(corrMapAvg)+0.2]);
            end

            % Plot correaltion heatmap vs duration of the signal taken into
            % consideration
            xlim([1 15]);xlabel('Time (minutes)'); ylabel('Correlation');
            plot(1:length(corrMapAvg),corrMapAvg,'k--','LineWidth',4); 
            patch([(1:length(corrMapAvg))' ;flipud((1:length(corrMapAvg))')],[(corrMapAvg-corrMapStd);  flipud((corrMapAvg+corrMapStd))],'r','FaceAlpha',0.2,'EdgeColor','none');
            xticks(1:15);title(figTitle); ylim([0.8 1.1]); yticks(0.8:0.05:1.1);
        end
        sgtitle(strrep([expDate ' Datafile: ' num2str(fileNum) ' - Correlation vs time '],'_','\_'));
%         f = gcf;
%         if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
%             [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
%         end
%         exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\DurationAnalysis_' bandName '_' num2str(fileNum) '.png'],'Resolution',300);
%         close gcf;
    end
end

%% Comparing 10 minute segment duration with 15 minute duration 
if strcmp(monkeyName,'CharlieSheen'); expDate = '08_30_2022'; iDate = 1;  datFileNum = datFileNumAll{iDate,1}; end 
if strcmp(monkeyName,'Whiskey'); expDate = '09_19_2022'; iDate = 2;  datFileNum = datFileNumAll{iDate,1}; end 

filterFlag = 0; binWidth = 60*fs;
for iFile = 1:length(datFileNum) % Get all datafiles for each experiment date
    fileNum = datFileNum(iFile);
    for iSet = 1:2
        clear probe pow freq figTitle chInCortex corrChAll tVals L chLen
        if (estChInCortexA{iDate}(fileNum,1)~= 0) && (estChInCortexB{iDate}(fileNum,1)~=0)
        switch iSet
            case 1
                probe = probe1{fileNum,iDate}(:,estChInCortexA{iDate}(fileNum,1):estChInCortexA{iDate}(fileNum,2));
                figTitle = 'Probe A';
                chInCortex = chInCortexProbeA{iDate};
                figure;
            case 2
                probe = probe2{fileNum,iDate};
                if size(probe,2)>1 % To separate the single channel probe
                    chVals = estChInCortexB{iDate}(fileNum,1):estChInCortexB{iDate}(fileNum,2);
                    chVals(ismember(chVals,badChProbeB)) = [];
                    probe = probe(:,chVals);
                end
                figTitle = 'Probe B';
                chInCortex = chInCortexProbeB{iDate};
        end
        else 
            continue;
        end
        if filterFlag
            probe = filtfilt(b,a,probe);
            bandName = 'Gamma';
        else
            bandName = 'WB'; % Wideband
        end

        % Get the 15 minute correlation heatmap
        corrChAll = corrcoef(probe,'rows','complete');
        tVals     = 10*60*fs; % Converting time to samples
        L         = size(probe,1); % Length of the signal collected
        chLen     = size(probe,2); % # of channels present in a recording
        if sum(tVals>L)>=1; tVals(tVals>L)=[]; end

        % Get the correlation heatmap for the recordings in a sliding
        % window with sliding window bin of 1 minute and variable bin
        % width
        clear corrMaps
        for iT = 1:length(tVals)
            w = tVals(iT);
            clear corrDat
            ind = 1;
            for iW = 1:binWidth:L-w+1
                if ~(size(probe,2)== 1)
                    corrDat(ind,:,:) = corrcoef(probe(iW:iW+w-1,:),'rows','complete');
                    corrMaps(:,ind) = corr(reshape(squeeze(corrDat(ind,:,:)),[chLen*chLen 1]),reshape(corrChAll,[chLen*chLen 1]));

                else % Single channel condition
                    corrDat(ind,:,:) = corrcoef(probe(iW:iW+w-1,:),'rows','complete');
                    corrMaps(:,ind) = 1;
                end
                ind = ind+1;
            end
        end

        corrMapAvg(fileNum,iSet) = mean(corrMaps,2); 
        corrMapStd(fileNum,iSet) = std(corrMaps);
        yLabel     = 'Correlation';
    end
end

% Plotting
corrMapAvgTmp = corrMapAvg; corrMapStdTmp = corrMapStd; 
corrMapAvgTmp([4 11],:) = []; corrMapStdTmp([4 11],:) = [];

figure;
for iSet = 1:2
    subplot(2,1,iSet); 
    plot(1:size(corrMapAvgTmp,1),corrMapAvgTmp(:,iSet),'*b','MarkerSize',10); hold on; 
    erBar = errorbar(1:size(corrMapAvgTmp,1),corrMapAvgTmp(:,iSet),corrMapStdTmp(:,iSet),corrMapStdTmp(:,iSet));
    xlabel('datafiles'); ylabel('Correlations');  xlim([ 1 12.5]);
    if iSet == 1; title('Probe A'); else; title('Probe B'); end 
    erBar.Color = [0 0 0]; 
    xticks((1:size(corrMapAvgTmp,1)));
    yline(1,'--');ylim([ 0.9 1.1]);
end 
sgtitle('Correlation between 10 minute heatmaps and 15 minute heatmaps ');
f = gcf;
exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\Corr10_15.png'],'Resolution',300);
close gcf;

%%  Anesthesia level analysis - 
alphaBand = [8 12];
betaBand = [15 30];
gammaBand = [30 90];

if strcmp(monkeyName,'Whiskey')
    expDate  = '10_17_2022';
    dateIdx  = strmatch(expDate,allDates);
    fileNum  = 5:8; % File Numbers containing data with varying anesthesia
    isoLevel = {'1.2', '2.25', '3.25', '1'}; % in %
    for iBand = 1:4
        if iBand == 4; figTitle = 'WB'; end
        clear b a corrProbes
        switch iBand
            case 1
                [b,a] = butter(3,alphaBand./(fs/2),'bandpass');
                figTitle = 'Alpha band';
            case 2
                [b,a] = butter(3,betaBand./(fs/2),'bandpass');
                figTitle = 'Beta band';
            case 3
                [b,a] = butter(3,gammaBand./(fs/2),'bandpass');
                figTitle = 'Gamma band';
        end
        for iC = 1: length(fileNum)
            x1 = probe1{fileNum(iC),dateIdx};
            y1 = probe2{fileNum(iC),dateIdx};
            y1(:,badChProbeB)= [];

            for iFile = 1:length(fileNum)
                chA = estChInCortexA{dateIdx}(fileNum(iFile));
                chB = estChInCortexB{dateIdx}(fileNum(iFile));
                if iBand ~=4
                    yA = probe1{fileNum(iFile),dateIdx};
                    yA(:,badChProbeB) = [];
                    xA = filtfilt(b,a,probe1{fileNum(iFile),dateIdx}(:,chA:chA+17));
                    yA = filtfilt(b,a,yA(:,chB:chB+17));
                else
                    xA = probe1{fileNum(iFile),dateIdx}(:,chA:chA+17);
                    yA = probe2{fileNum(iFile),dateIdx}(:,chB:chB+17);
                end
                corrProbesMean(iFile,:) = mean(corr(xA,yA),2,'omitnan'); % Pairwise correlations between probe A and probe B

                % Correlations for the two probes under different levels of
                % anesthesia
                x = probe1{fileNum(iFile),dateIdx};
                y = probe2{fileNum(iFile),dateIdx};  y(:,badChProbeB) = [];

                if iBand ~=4
                    xFilt = filtfilt(b,a,x);
                    yFilt = filtfilt(b,a,y);
                    xFirst = filtfilt(b,a,x1);
                    yFirst = filtfilt(b,a,y1);
                else
                    xFilt  = x;  yFilt   = y;
                    xFirst = x1; yFirst = y1;

                end

                corrProbeA(iBand,iC,iFile,:) = corr2(imgaussfilt(corr(xFirst),1), imgaussfilt(corr(xFilt),1));
                corrProbeB(iBand,iC,iFile,:) = corr2(imgaussfilt(corr(yFirst),1), imgaussfilt(corr(yFilt),1));
            end
        end

        figure; boxplot(corrProbesMean','Notch','off','Labels',isoLevel); ylim([ -0.1 1]); hold on;
        lines = findobj(gcf,'type','line','Tag','Median');
        xVal = mean(vertcat(lines.XData),2);
        yVal = vertcat(lines.YData);
        plot(xVal,yVal(:,1),'k--','LineWidth',1);
        ylabel('Correlation values'); xlabel(' Isoflurane levels (%)');
        title(['Pairwise correlations at different levels of isoflurane - ' figTitle]);

%         figure; plot(corrProbeA,'r*--','MarkerSize',8,'LineWidth',1); hold on; 
%         plot(corrProbeB,'b*--','MarkerSize',8,'LineWidth',1); legend('Probe A','ProbeB','Location','northeast');
%         title(['Intra-probe correlations at different levels of isoflurane - ' figTitle]);
%         xticks(1:4); xticklabels(isoLevel); ylim([ -0.1 1.1]); yticks(-0.1:0.1:1.1);
%         xlabel(' Isoflurane levels (%)'); ylabel('Correlation between intra probe heatmaps'); 

    end

    colorVals = hsv(length(fileNum));
    for iBand = 1:4
        switch iBand
            case 1
                figTitle = 'Alpha band';
            case 2
                figTitle = 'Beta band';
            case 3
                figTitle = 'Gamma band';
            case 4
                figTitle = 'Wideband';
        end
        for iSet = 1:2
            switch iSet
                case 1
                    probeVals = corrProbeA;
                    probeTitle = 'Probe A';
                case 2
                    probeVals = corrProbeB;
                    probeTitle = 'Probe B';
            end
            figure; 
            for iC = 1:length(fileNum)
               plot(squeeze(probeVals(iBand,iC,:)),'*--','Color',colorVals(iC,:,:),'MarkerSize',8,'LineWidth',1); hold on;
            end
            xticks(1:4); xticklabels(isoLevel); ylim([ -0.2 1.2]); xlabel('Anesthesia levels (%)'); ylabel('Correlation values');
            title(['Intraprobe correlations for ' probeTitle ' - ' figTitle]); 
            legend(isoLevel,'Location','southeast');
        end
    end
end


 

%% OLD CODE %% 
%  Find out which channels  are inside the cortex - get the intrachannel correlation for probe 1 and probe 2 - Under progress...
clear b a
plotCorrMapFlag = 0;
plotRelCorrFlag = 0;
plotAvgCorrFlag = 0;
plotRMSFlag     = 0;
plotPowFlag     = 1;
filterFlag      = 0;
fs              = 1e3;
alphaBand       = [8 12];
betaBand        = [15 30];
gammaBand       = [30 90];
[b,a]           = butter(3,[30 90]./(fs/2),'bandpass');

for iDate = 2%1:size(allDates,1)
    datFileNum = datFileNumAll{iDate,1};
    expDate = allDates(iDate,:);
    for iFile = 10%1:length(datFileNum)
        for iSet = 1:2
            clear probe pow freq figTitle chInCortex
            switch iSet
                case 1

                    probe = probe1{iFile,iDate};
                    figTitle = 'Probe A';
                    chInCortex = chInCortexProbeA{iDate};
                case 2
                    probe = probe2{iFile,iDate};
                    if size(probe,2)>1 % To separate the single channel probe
                        probe(:,badChProbeB)= [];
                    end
                    figTitle = 'Probe B';
                    chInCortex = chInCortexProbeB{iDate};
            end

            if filterFlag % Filtering gamma band
                probe = filtfilt(b,a,probe);
                bandName = 'Gamma';
            else
                bandName = 'WB'; % Wideband
            end

            if plotCorrMapFlag && chInCortex(iFile)~=1  % Plot intra channel correlation heatmaps/correlogram
                corrChannels = corrcoef(probe,'rows','complete');
                figure; imagesc(imgaussfilt(corrChannels,1)); colormap jet; axis image; colorbar; hold on;
                caxis([0 1]); xticks(1:length(channels)); yticks(1:length(channels));
                sgtitle(strrep([datFileName(1:8) ' ' num2str(iFile) ' ' figTitle],'_','\_'));

                f = gcf;
                if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
                    [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
                end
                exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\' bandName '_heatmap_' num2str(iFile) '_' num2str(iSet) '.png'],'Resolution',300);
                close gcf;
            end
            % Plot the relative correlation
            if plotRelCorrFlag == 1
                corrChannels = corrcoef(probe,'rows','complete');
                corrColWise = mean(corrChannels,1);
                subplot(1,2,iSet);
                barh(corrColWise);
                yticks(1:32); xlabel('Relative correlation'); ylabel('Channels');
                title(figTitle);
                figSaveName = 'RelCorr';
            end

            %         chOutlier = find(isoutlier(rms));
            %         [~,chIdx] = min(abs(chInCortex(iFile)-chOutlier));
            %         if ~isempty(newChInCortex)
            %             newChInCortex = chOutlier(chIdx)+1;
            %         else
            %             newChInCortex = chInCortex(iFile);
            %         end

            % Get the ratio of RMS values and adjust the channels that are inside vs outisde the cortex
            if plotRMSFlag && chInCortex(iFile)~=1
                if chInCortex(iFile)~=1
                    probe = [probe(:,1) diff(probe')'];
                else
                    probe = probe;
                end

                rms = sqrt(mean(probe.^2))./sqrt(mean(probe(:,chInCortex(iFile):size(probe,2)).^2,[1 2]));
                upperBound = mean(rms)+1*std(rms);
                newChInCortex = find(rms>upperBound, 1);
                subplot(2,1,iSet);
                plot(1:size(probe,2),rms,'r*','MarkerSize',8); hold on; yline(1);
                if iSet == 1; ylabel('RMS/(RMS of channels in cortex)'); end
                xlabel('Channels'); xticks(1:size(probe,2));
                title(['Datafile: ' num2str(iFile) ' ' figTitle ]); ylim([ 0 10]);
                figSaveName = 'RMS';
            end

            % Replace the channels that are in cortex if already not changed
            %         switch iSet
            %             case 1
            %                 if chInCortexProbeA(iFile) ~= newChInCortex
            %                     chInCortexProbeA(iFile) = newChInCortex;
            %                 end
            %             case 2
            %                 if chInCortexProbeB(iFile) ~= newChInCortex
            %                     chInCortexProbeB(iFile) = newChInCortex;
            %                 end
            %         end

            % Plotting the Power spectrum
            if plotPowFlag
                if chInCortex(iFile)~=1
                    probe = [probe(:,1) diff(probe')'];
                else
                    probe = probe;
                end
                figure;
                for iElec = 1:size(probe,2)
                    [freq,pow(:,iElec)] = getFFT(probe(:,iElec),fs,0);
                    %             plot(freq,10.*log10(pow(:,iElec).^2));
                    if chInCortex(iFile)~=1
                        subplot(4,8,iElec);
                    end
                    plot(freq,pow(:,iElec)); xlim([0 50]); hold on;
                    xlabel('Frequency');ylabel('Power');
                    title(['Channel: ' num2str(iElec)]);
                    ylim([0 20]);
                end

                sgtitle(['Datafile: ' num2str(iFile) ' ' figTitle ]);

                % Doing PCA  to determine the bounds...
                %             trainVals = 3:size(probe,2)-2;
                %             testVals = [1 2 size(probe,2)-1 size(probe,2)];

                % Centering the data
                centeredPow =  pow - mean(pow,2);
                %             centeredTest =   bsxfun(@minus, pow(:,testVals),mean(pow(:,testVals),2));
                %
                %             [pcVecs,pcScores,eigVals] = pca(centeredPow(1:9023,:)); %PCA
                %             percentVar = (eigVals./sum(eigVals)).*100;


                %             projData = centeredPow(1:9023,:)*pcVecs(:,1:2);
                %             reconTrainData = projTrainData*pcVecs(:,1:7)';

                %             projTestData = centeredTest*pcVecs(:,1:7); %Test data in projected PC space
                %             reconTestData = projTestData*pcVecs(:,1:7)';

                %             r2Train = compute_r2(centeredTrainTarg,reconTrainData);
                %             r2Test = compute_r2(centeredTestTarg,reconTestData);
                %             figure;
                %             for iElec = 1:size(probe,2)
                %                 subplot(8,4,iElec);
                % %                 plot(projData(:,iElec)+ (iElec-1)*500);
                % %                 plot(meanElecData{1,eventID}(:,iElec)+ ((iElec-1)*100));
                % %                 plot(P(:,iElec)+ (iElec-1)*100);
                %                   plot(projData(:,iElec)); ylim([ -10 30]);
                %                 hold on;
                %             end

            end

        end

        %         if plotRMSFlag || plotRelCorrFlag
        %             sgtitle(strrep([datFileName(1:8) ' ' num2str(iFile) ],'_','\_'));
        %             f = gcf;
        %             if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
        %                 [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
        %             end
        %             exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\' figSaveName '_' bandName '_' num2str(iFile) '.png'],'Resolution',300);
        %             close gcf;
        %         end
    end
end



%% Looking at time-wise correlation heatmap
% Get a set of 15 heat maps - 1 minute,2 minute, 3 minute.... 15 minutes
binWidth = 60*fs;  % 1 minute bin width
for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        for iSet = 1:2
            clear probe pow freq figTitle chInCortex corrChAll tVals L chLen
            switch iSet
                case 1
                    probe = probe1{iFile,iDate};
                    figTitle = 'Probe A';
                    chInCortex = chInCortexProbeA{iDate};

                case 2
                    probe = probe2{iFile,iDate};
                    if size(probe,2)>1 % To separate the single channel probe
                        probe(:,badChProbeB)= [];
                    end
                    figTitle = 'Probe B';
                    chInCortex = chInCortexProbeB{iDate};
            end

            % Get the 15 minute correlation heatmap
            corrChAll = corrcoef(probe,'rows','complete');
            tVals = (1:15).*60*fs; % Converting to samples
            L = size(probe,1); % Length of the signal collected
            chLen = size(probe,2);
            % Get the correlation heatmap for the recordings in a sliding
            % window with sliding window bin of 1 minute and variable bin
            % width
            figure;
            if sum(tVals>L)>=1; tVals(tVals>L)=[]; end
            clear corrMaps winAvgCorr
            for iT = 1:length(tVals)
                w = tVals(iT);
                clear corrDat
                ind = 1;
                for iW = 1:binWidth:L-w+1
                    if ~(size(probe,2)== 1)
                        corrDat(ind,:,:) = corrcoef(probe(iW:iW+w-1,:),'rows','complete');
                        corrMaps{iT,1}(:,ind) = corr(reshape(squeeze(corrDat(ind,:,:)),[chLen*chLen 1]),reshape(corrChAll,[chLen*chLen 1]));

                    else % Single channel condition
                        corrDat(ind,:,:) = corrcoef(probe(iW:iW+w-1,:),'rows','complete');
                        corrMaps{iT,1}(:,ind) = 1;
                    end
                    ind = ind+1;
                end
                timeWiseCorr{iT,1} = squeeze(median(corrDat,1,'omitnan'));
                subplot(3,5,iT);
                imagesc(imgaussfilt(timeWiseCorr{iT,1},1)); colormap jet; caxis([0 1]);
                xticks(1:5:chLen); yticks(1:5:chLen);
                title(['Time - ' num2str(w/(60*fs)) ' min']);
            end
            sgtitle(strrep([expDate ' Datafile' num2str(iFile) ' ' figTitle],'_','\_'));
            f = gcf;
            if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
                [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
            end
            exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\TimeWiseHeatMap' num2str(iFile) '_' figTitle '.png'],'Resolution',300);
            close gcf;
        end
    end
end

%% Cumulative heatmap from 1 minute to 15 minutes
filterFlag = 1; saveFlag = 0;
gammaBand = [30 90];
[b,a] = butter(3,gammaBand./(fs/2),'bandpass');

for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        for iSet = 1:2
            clear probe pow freq figTitle chInCortex corrChAll tVals L chLen
            switch iSet
                case 1
                    if ~saveFlag;figure;end
                    probe = probe1{iFile,iDate};
                    figTitle = 'Probe A';
                    chInCortex = chInCortexProbeA{iDate};

                case 2
                    probe = probe2{iFile,iDate};
                    if size(probe,2)>1 % To separate the single channel probe
                        probe(:,badChProbeB)= [];
                    end
                    figTitle = 'Probe B';
                    chInCortex = chInCortexProbeB{iDate};
            end

            if filterFlag % Filtering gamma band
                probe = filtfilt(b,a,probe);
                bandName = 'Gamma';
            else
                bandName = 'WB'; % Wideband
            end

            % Get the 15 minute correlation heatmap
            corrChAll = corrcoef(probe,'rows','complete');
            corrChAllR = reshape(corrChAll,[size(corrChAll,1)*size(corrChAll,2) 1]);
            tVals = (1:15).*60*fs; % Converting to samples
            L = size(probe,1); % Length of the signal collected
            chLen = size(probe,2);
            % Get the correlation heatmap for the recordings in a sliding
            % window with sliding window bin of 1 minute and variable bin
            % width

            if sum(tVals>L)>=1; tVals(tVals>L)=[]; end
            clear corrMaps winAvgCorr corrDat compareCorr
            for iT = 1:length(tVals)
                w = tVals(iT);
                %                 clear corrDat

                if ~(size(probe,2)== 1)
                    corrDat(iT,:,:) = corrcoef(probe(1:w,:),'rows','complete');
                    compareCorr(iT,:) = corr(reshape(squeeze(corrDat(iT,:,:)),[size(corrDat,2)*size(corrDat,3) 1]),corrChAllR);
                    if saveFlag
                        subplot(3,5,iT);
                        imagesc(imgaussfilt(squeeze(corrDat(iT,:,:)),1)); colormap jet; caxis([0 1]);
                        xticks(1:5:chLen); yticks(1:5:chLen);
                        title(['Time - ' num2str(w/(60*fs)) ' min']);
                    end
                else % Single channel condition
                    corrDat(iT,:,:) = corrcoef(probe(1:w,:),'rows','complete');
                    compareCorr(iT,:) = 1;
                    continue;
                end
            end
            %             avgTime = squeeze(mean(corrDat,[2 3],'omitnan'));
            subplot(2,1,iSet); plot(1:length(compareCorr),compareCorr,'k--o','MarkerSize',5); hold on;
            xlabel('Time'); ylabel('Correlation');title([figTitle ' - ' bandName]);
            xlim([ 1 length(compareCorr)]); xticks(1:2:length(compareCorr));
            ylim([ min(compareCorr)-0.1 max(compareCorr)+0.1]);

            if saveFlag
                if ~(size(probe,2)== 1)
                    sgtitle(strrep([expDate ' Datafile' num2str(iFile) ' ' figTitle],'_','\_'));
                    f = gcf;
                    if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
                        [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
                    end
                    exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\CumulativeTimeHeatMap_' num2str(iFile) '_' bandName '_' figTitle '.png'],'Resolution',300);
                    close gcf;
                end
            end
        end
        sgtitle(strrep([expDate ' Datafile' num2str(iFile)],'_','\_'));
        f = gcf;
        if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'])
            [~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results']);
        end
        exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\CorrVsCumTime_' num2str(iFile) '_' bandName '.png'],'Resolution',300);
        close gcf;

    end
end

%% Pairwise correlation heatmaps between probe A and probe B
filterFlag = 1; saveFlag = 0;
% gammaBand = [30 90];
alphaBand       = [8 12];
[b,a] = butter(3,alphaBand./(fs/2),'bandpass');

binWidth   = 60*fs;
for iDate = 1:size(allDates,1)
    clear expDate datFileNum
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    if strcmp(expDate,'08_08_2022'); continue; end

    for iFile = 1:length(datFileNum)
        if strcmp(expDate,'09_19_2022') && iFile == 4; continue; end
        probeA = probe1{iFile,iDate};
        probeB = probe2{iFile,iDate};

        if filterFlag
            probeA = filtfilt(b,a,probeA);
            probeB = filtfilt(b,a,probeB);
            bandName = 'Alpha';
        else
            bandName = 'WB';
        end

        probeB(:,badChProbeB)= [];
        corrPairwise = corr(probeA,probeB);
        tVals    = (1:15).*60*fs; % Converting time to samples
        L        = size(probeA,1); % Length of the signal collected
        if sum(tVals>L)>=1; tVals(tVals>L)=[]; end

        % Get the correlation heatmap for the recordings in a sliding
        % window with sliding window bin of 1 minute and variable bin
        % width
        clear corrMaps winAvgCorr
        figure;
        for iT = 1:length(tVals)
            w = tVals(iT);
            clear corrDat
            ind = 1;
            for iW = 1:binWidth:L-w+1
                corrDat(ind,:,:) = corr(probeA(iW:iW+w-1,:),probeB(iW:iW+w-1,:));
                corrMaps{iT,1}(:,ind) = corr(reshape(squeeze(corrDat(ind,:,:)),[size(corrDat,2)*size(corrDat,3) 1]), ...
                    reshape(corrPairwise,[size(corrPairwise,1)*size(corrPairwise,2) 1]));
                ind = ind+1;
            end
            winAvgCorr{iT,1} = mean(corrDat,[2,3],'omitnan');
        end
        corrMapAvg = cellfun(@mean,winAvgCorr);
        corrMapStd = cellfun(@std,winAvgCorr);

        for iL = 1:length(tVals)
            plot(iL.*ones(size(winAvgCorr{iL,1})),winAvgCorr{iL,1},'b.','MarkerSize',8); hold on;
        end
        plot(1:length(corrMapAvg),corrMapAvg,'k--','LineWidth',2);
        patch([(1:length(corrMapAvg))' ;flipud((1:length(corrMapAvg))')],[(corrMapAvg-corrMapStd);  flipud((corrMapAvg+corrMapStd))],'r','FaceAlpha',0.2,'EdgeColor','none');
        xticks(1:2:15);
        xlabel('Time'); ylabel('Average pairwise correlation');title(strrep([expDate ' Datafile ' num2str(iFile) ' - ' bandName],'_','\_'));
        ylim([ min(corrMapAvg)-0.2 max(corrMapAvg)+0.2]);

        f = gcf;
        exportgraphics(f,['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\AvgPaircorrVsTime_' num2str(iFile) '_' bandName '.png'],'Resolution',300);
        close gcf;
    end
end

%% Transition -inside vs outside - comparing one channel inside cortex with that of others and keeping controls


%%  OLD CODE %%
% To find out which probes are inside the cortex - to do the following-
% 1. Get the power spectrum for all the chaneks
% 2. Get the RMS value of each channel and divide it over the RMS value of
% channels inside the cortex to get the RMS ratio
% 3. Get the SNR with the noise being expected value of signals inside the
% cortex (or channels that are inside the cortex)
for iFile = 5%:length(datFileNum)
    for iSet = 1:size(uniqueCh,1)
        clear probe pow freq figTitle chInCortex
        switch iSet
            case 1
                probe = probe1{iFile};
                figTitle = 'Probe A';
                chInCortex = chInCortexProbeA:size(probe,2);
            case 2
                probe = probe2{iFile};
                probe(:,badChProbeB)= [];
                figTitle = 'Probe B';
                chInCortex = chInCortexProbeB:size(probe,2);
        end

        %         % Plotting the Power spectrum
        %         ax = figure;
        %         for iElec = 1:size(probe,2)
        %             [freq,pow(:,iElec)] = getFFT(probe(:,iElec),fs,0);
        %             ax(iElec) = subplot(4,8,iElec);
        % %             plot(freq,10.*log10(pow(:,iElec).^2));
        %             plot(freq,pow(:,iElec)); xlim([0 70]);
        %             xlabel('Frequency');ylabel('Power');
        %             title(['Channel: ' num2str(iElec)]);
        %             ylim([0 20]);
        %         end
        %         sgtitle(['Datafile: ' num2str(iFile) ' ' figTitle ]);
        %         linkaxes(ax,'y');

        % Plotting the RMS value
        rms = sqrt(mean(probe.^2))./sqrt(mean(probe(:,chInCortex).^2,[1 2]));
        figure;
        plot(1:size(probe,2),rms,'r*','MarkerSize',15);
        ylabel('RMS values with respect to the RMS of channels in cortex'); xlabel('Channels'); xticks(1:size(probe,2));
        title(['Datafile: ' num2str(iFile) ' ' figTitle ]);

        %         % Plotting the SNR, given we know the channels that are within the
        %         % cortex.
        %         snrVals = 10.*log10(sqrt(mean(probe.^2))./ sqrt(mean(probe(:,chInCortex).^2,[1 2])));
        %         figure;
        %         plot(1:size(probe,2),snrVals,'k*','MarkerSize',15);
        %         ylabel('SNR (in dB)'); xlabel('Channels');
        %         title(['Datafile: ' num2str(iFile) ' ' figTitle ]); xticks(1:32);
    end
end

%% Get the powers separately based on channels that are inside and outside the cortex
badChProbeB = [14 22];
chInCortexProbeA = [10 5 8 9 5 11 11 4 4 8 8 8];
chInCortexProbeB = [10 10 10 10 10 8 8 8 8 8 9 5];
for iFile = 1:length(datFileNum)
    for iSet = 1:size(uniqueCh,1)
        clear probe pow freq figTitle chInCortex colorVal
        switch iSet
            case 1
                probe = probe1{iFile};
                figTitle = 'Probe A';
                chInCortex = chInCortexProbeA(iFile):size(probe,2);
            case 2
                probe = probe2{iFile};
                probe(:,badChProbeB)= [];
                figTitle = 'Probe B';
                chInCortex = chInCortexProbeB(iFile):size(probe,2);
        end

        % Plotting the Power spectrum
        figure;
        for iElec = 1:size(probe,2)
            [freq,pow(:,iElec)] = getFFT(probe(:,iElec),fs,0);
            if iElec< chInCortex(1)
                colorVal = [ 1 0 0 0.2];
            else
                colorVal = [0 0 1 0.2];
            end
            %             plot(freq,10.*log10(pow(:,iElec).^2));
            plot(freq,10.*log10(pow(:,iElec).^2),'Color',colorVal); xlim([0 50]); hold on;
            xlabel('Frequency');ylabel('Power (dB)');
            %             title(['Channel: ' num2str(iElec)]);
            ylim([-100 60]);
        end
        sgtitle(['Datafile: ' num2str(iFile) ' ' figTitle ]);

    end
end

%% Plotting pairwise correlation for all the recorded sites
for iFile = 1:length(datFileNum)
    for iSet = 1:size(uniqueCh,1)
        clear probe pow freq figTitle chInCortex colorVal
        switch iSet
            case 1
                probe = probe1{iFile};
                figTitle = 'Probe A';
                chInCortex = chInCortexProbeA(iFile):size(probe,2);
                probe1New{iFile,1} = probe(:,chInCortex);
            case 2
                probe = probe2{iFile};
                probe(:,badChProbeB)= [];
                figTitle = 'Probe B';
                chInCortex = chInCortexProbeB(iFile):size(probe,2);
                probe2New{iFile,1} = probe(:,chInCortex);
        end
    end
end
%% Mean correlation wrt distance between sites
for iFile  = 1:length(datFileNum)
    clear xA yA
    xA = probe1New{iFile};
    yA = probe2New{iFile};
    %     if size(x,2)>size(y,2)
    yA = [yA NaN(size(yA,1), 32-size(yA,2))];
    %     else
    xA = [xA NaN(size(xA,1), 32-size(xA,2))];
    %     end
    meanCorrProbes(:,iFile) = mean(corr(xA,yA),[1 2],'omitnan');
    corrProbes(:,:,iFile) = corr(xA,yA);
end

connSites = logical([1 1 1 0 0 1 0 1 0 0 1 1]);
distSitesAll = [3.6 3 0.8 4.13 5.25 3.77 3.8 3.2 5 2.6 4.36 4.38];
distConnSites = distSitesAll(connSites);
distNCSites = distSitesAll(~connSites);
meanCorrConnSites = meanCorrProbes(:,connSites);
meanCorrNCSites = meanCorrProbes(:,~connSites);
figure; plot(distConnSites,meanCorrConnSites,'r*'); hold on;
plot(distNCSites,meanCorrNCSites,'bo'); xlabel('Distance (mm)'); ylabel('Mean pairwise correlations'); legend({'Connected'; ...
    'Unconnected'});
figure; boxplot([meanCorrConnSites' [meanCorrNCSites'; zeros(length(meanCorrConnSites)-length(meanCorrNCSites),1)]],'Notch','off','Labels',{'Connected', 'Not connected'}) ;%hold on; boxplot(meanCorrNCSites);
ylabel('Correlation values');

%%

% Comparing the correlations between probe 1 and probe 2 for different datafiles
%  figure; boxplot(squeeze(corrProbes(:,1,1:4));%,'Notch','off','Labels',{'Datafile0001', 'Datafile0002','Datafile0003','Datafile0004'}); hold on;
%  plot(repmat(1:4, [32 1]),squeeze(corrProbes(:,1,1:4)),'g.','Markersize',10);
ylabel('Correlation between probe 1 and probe 2');


%% Checking which channels are good or bad
for iFile = 1:length(datFileNum)
    for iSet =2% 1:size(uniqueCh,1)
        clear probe pow freq figTitle chInCortex colorVal
        switch iSet
            case 1
                probe = probe1{iFile};
                figTitle = 'Probe A';
                chInCortex = chInCortexProbeA(iFile):size(probe,2);
            case 2
                probe = probe2{iFile};
                probe(:,badChProbeB)= [];
                figTitle = 'Probe B';
                chInCortex = chInCortexProbeB(iFile):size(probe,2);
        end
        for iCh = 1:size(probe,2)
            c14(iCh,1) = corr(probe(:,14),probe(:,iCh));
            c22(iCh,1) = corr(probe(:,22),probe(:,iCh));
        end

    end
end



%%
%             if strcmp(lfpLabel{iElec}(end-2:end),'lfp')
%             else
%                 if strcmp(lfpLabel{iElec},'analog 2') % analog 2 has the single probe
%                     clear probe2Temp
%                     [~, ~, probe2Temp] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
%                     [~, analogInfo2] = ns_GetAnalogInfo(hFile, lfpEntityID);
%                     fs_ns5  = analogInfo2.SampleRate;
%                     % Downsample the analog 2 channel from 30kHz to 1kHz and
%                     % bandpass filter across 1- 250 Hz
%                     probe2{iFile} = downsample(probe2Temp,fs_ns5/fs);
%                     probe2{iFile} = filtfilt(b,a,probe2{iFile});
%
% %                 else % analog 1 has EEG
% %                     [~, ~, eeg{iFile}] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
%                 end
%             end



%         end
%
%
%
% %         elecID  = unique([hFile.Entity(lfpList).ElectrodeID]); % Get the unique lfp channels
% %         elecIDTemp = elecID; % Keeping the original copy of the electrode ID
% %         analogElecs = find(elecID>max(lfpList)); % To get the analog channels containing the EEG and/or single probe
% %         if ~isempty(analogElecs); elecID(analogElecs) = analogElecs + 32; end
%
%         % Obtaining the entity labels and organizing according to the
%         % channel map
% %         lfpLabel = {entityInfo(lfpList(elecID(1:end-1))).EntityLabel};
%         channels = length(lfpLabel);
%         chNum = cell2mat(cellfun(@(x) x(2:3),lfpLabel','un',0));
%         chNum = str2num(chNum(1:end-2,:));
%         [~,elecIDTemp(1:end-2)] = sort(chNum);  % Rearranging the channels based on the channel map
%
%         %% Get the LFP for all the channels in sorted order
%         for iElec = 1:channels
%             clear elecEntityID lfpEntityID lfpCount
%             elecEntityID  = find([hFile.Entity(:).ElectrodeID] == elecIDTemp(iElec));
%             lfpEntityID   = elecEntityID(1);
%             lfpCount      = entityInfo(lfpEntityID).ItemCount;
%
%             if ~exist('analogInfo','var')
%                 [~, analogInfo] = ns_GetAnalogInfo(hFile, lfpEntityID);
%                 timeStamps = hFile.FileInfo(hFile.Entity(lfpEntityID).FileType).TimeStamps;
%                 fs  = analogInfo.SampleRate;
%                 analogTimes   = (0:(sum(timeStamps(:,end)))-1)' ./ fs;
%                 clear b a
%                 [b,a] = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering  across 1-250 Hz
%             end
%
%             if strcmp(lfpLabel{iElec}(end-2:end),'lfp')
%                 [~, ~, probe1{iFile}(:,iElec)] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount); % Get LFP data
%                  probe1{iFile} =    filtfilt(b,a,probe1{iFile});
%             else
%                 if strcmp(lfpLabel{iElec},'analog 2') % analog 2 has the single probe
%                     clear probe2Temp
%                     [~, ~, probe2Temp] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
%                     [~, analogInfo2] = ns_GetAnalogInfo(hFile, lfpEntityID);
%                     fs_ns5  = analogInfo2.SampleRate;
%                     % Downsample the analog 2 channel from 30kHz to 1kHz and
%                     % bandpass filter across 1- 250 Hz
%                     probe2{iFile} = downsample(probe2Temp,fs_ns5/fs);
%                     probe2{iFile} = filtfilt(b,a,probe2{iFile});
%
% %                 else % analog 1 has EEG
% %                     [~, ~, eeg{iFile}] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
%                 end
%             end
%
%         end
%     end




%% Isolate the channels that are outside the cortex and take the channels that are inside the cortex only
for iFile = 1:length(datFileNum)
    for iSet = 1:size(uniqueCh,1)
        clear probe
        switch iSet
            case 1
                probe = probe1{iFile,1};
            case 2
                probe = probe2{iFile,1};
        end
        figure;
        corrChannels(:,:,iFile) = corrcoef(probe);
        imagesc(imgaussfilt(corrChannels(:,:,iFile),1)); colormap jet; axis image; colorbar
        caxis([ min(corrChannels(:,:,iFile),[],'all') max(corrChannels(:,:,iFile),[],'all')])
        xticks(1:length(channels));
        yticks(1:length(channels));
        title(strrep([datFileName(1:8) ' ' num2str(iFile)],'_','\_')); drawnow;

        chInput(1) = input('Enter the channel that is the most superficial ');
        chInput(2) = input('Enter the channel that is the deepest ');
        probe = probe(:,chInput(1):chInput(2));
        otherChFlag = input('Do you want to remove other channels? 1- Yes, 0 - No ');
        if otherChFlag
            otherCh = input('Which other channel(s) do you want to remove (put the channels in square brackets) ');
            otherCh = otherCh - chInput(1);
            probe = probe(:,[1:otherCh(1)-1 otherCh(end)+1:end]);
        end
        switch iSet
            case 1
                probe1{iFile,1} = probe;
            case 2
                probe2{iFile,1} = probe;
        end
        close gcf;
    end
end


% Get the average pairwise correlation values between probe 1 and probe 2
for iFile  = 1:length(datFileNum)
    clear xA yA
    xA = probe1{iFile};
    yA = probe2{iFile};
    if size(xA,2)>size(yA,2)
        yA = [yA NaN(size(yA,1), size(xA,2)-size(yA,2))];
    else
        xA = [xA NaN(size(xA,1), size(yA,2)-size(xA,2))];
    end
    meanCorrProbes(:,iFile) = mean(corr(xA,yA),[1 2],'omitnan');
end
d = [4 2 1 6 7 5 6 4 8 5]; % Rough estimate - in mm  - to verify thissss



%%
% Get the correlation values between the 32 probes vs the single probe for the resting state
figure;
for iFile = 1:length(datFileNum)
    %     subplot(2,2,iFile);
    corrProbes(:,:,iFile) = corr(probe1{iFile},probe2{iFile});
    %     plot(corrProbes(:,iFile));%,'.k','MarkerSize',10);
    %     xlabel('Channels in probe 1');
    %     ylabel('Correlation values');
    %     ylim([-0.25 0.25]);
    %     xlim([ 0 33]); xticks(1:33);
    %     title(strrep([datFileName num2str(iFile)],'_','\_'));

end

%% Comparing the correlations between probe 1 and probe 2 for different datafiles
figure; boxplot(squeeze(corrProbes(:,30,1:4)),'Notch','off','Labels',{'Datafile0001', 'Datafile0002','Datafile0003','Datafile0004'}); hold on;
%  plot(repmat(1:4, [32 1]),squeeze(corrProbes(:,1,1:4)),'g.','Markersize',10);
ylabel('Correlation between probe 1 and probe 2');


%% Check for artifact/noise in the LFP signals by calculating the signal to noise ratio of the LFP
% We are taking the most superficial channel (channel 1) and considering it
% as noise and calculating the SNR
% For datafile0002
iFile = 2;
snr = 10.*log10(sqrt(mean(probe1{2}.^2))./ sqrt(mean(probe1{2}(:,1).^2)));
figure; plot(snr,'k*','MarkerSize',13);  hold on;
m = mean(snr(11:32)); s = std(snr(11:32)); ylabel('SNR'); xlabel('Channels');
yline(m-s,'b-','Mean+std');
yline(m-0.5*s,'r-','Mean-0.5std');
xline(11,'-','Channels inside cortex','FontSize',13);
annotation('textarrow',[0.4 0.5],[ 0.85 0.85],'LineWidth',1);


%% Filter out the LFP based on different frequency bands

