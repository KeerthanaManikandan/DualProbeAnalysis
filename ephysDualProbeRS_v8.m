%% ephysDualProbeRS_v8
% This function performs analysis on LFP recorded simultaneously from two linear electrode arrays.
% January 14, 2026 - Keerthana Manikandan
% This code performs the following for ONE MONKEY (for combined analyis,
% check this script- combinedAnalysisEphysDualProbeRS.m)
% Check ephysDualProbeRS_v7 for previous iterations

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
rmpath(genpath([commonDir '\Codes\ISOI_Ephys']));
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
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate,patchInfo,pairClass,corticalAreaProbeA,...
    corticalAreaProbeB] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']);

% Get the connectivity and distance between pairs
clear distSites connSites greenMapRef
[distSites,connVals,refSites,movSites,greenMapRef] = ...
    getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);


distSitesAll = NaN(size(goodRuns)); connValsAll = NaN(size(goodRuns));
heartRateValsAll = NaN(size(goodRuns)); anesthesiaValsAll = NaN(size(goodRuns));

siteProbeA = cell(size(goodRuns)); siteProbeB = cell(size(goodRuns));


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
    datFileNameAll, datFileNumAll,serverPath,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);

%% Save bad channels, bad times in LFP mat file
% for iDate = 1: size(allDates,1)
%     clear expDate
%     expDate = allDates(iDate,:);
%     saveFolder = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
% 
%     for iRun = 1:length(datFileNumAll{iDate,1})
%         clear fileNum zL pL kL
%         fileNum = datFileNumAll{iDate,1}(iRun);  
% 
%         % Get the name of stored file
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
%         badElecA = badElecAall{fileNum,iDate};
%         badElecB = badElecBall{fileNum,iDate};
%         chCortexA = estChInCortexA{iDate}(fileNum,:);
%         chCortexB = estChInCortexB{iDate}(fileNum,:);
%         badTimes = allBadTimes{fileNum,iDate};
% 
%         save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'badElecA','badElecB',...
%             'chCortexA','chCortexB','badTimes','-append');
%     end
% end

%% Get all pairwise correlations
bandLabels = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking'};
timeLabels = {'Time series','Power','Infraslow'};
[allVars] = getDualProbeCorrelations(monkeyName, hemisphere, allDates, datFileNumAll,allProbeData,allBadTimes,...
    badElecA,badElecB,estChInCortexA,estChInCortexB,connValsAll,distSitesAll,goodRuns);

%% Get anesthesia recordings

dateVal = 3; runVal = 5:8; 
withinA = NaN(31,31,4,'single');
withinB = NaN(30,30,4,'single');

for iRun = 1:4
    
    chA = estChInCortexA{dateVal}(runVal(iRun),:);
    chB = estChInCortexB{dateVal}(runVal(iRun),:) ;

    probeA = allProbeData{runVal(iRun),dateVal}.probe1Ch; 
    probeB = allProbeData{runVal(iRun),dateVal}.probe2Ch; 


    probeA(:,badElecA{runVal(iRun),dateVal}) = []; 
    probeB(:,badElecB{runVal(iRun),dateVal}) = []; 

    % probeA = probeA(:,chA(1):chA(2));
    % probeB = probeB(:,chB(1):chB(2)); 

    probeA(allBadTimes{runVal(iRun),dateVal},:) = []; 
    probeB(allBadTimes{runVal(iRun),dateVal},:) = []; 

    % Get between and within probe correlations

    withinA(:,:,iRun) = corr(probeA,'Rows','complete');
    withinB(1:30,1:30,iRun) = corr(probeB,'Rows','complete'); 

    btwAB(:,:,iRun)   = corr(probeA,probeB,'Rows','complete');   

    meanCorrAB(iRun,1) = median(btwAB(chA(1):chA(2),chB(1):chB(2),iRun),[1 2],'omitnan');
    meanCorrA(iRun) = median(withinA(chA(1):chA(2),chA(1):chA(2),iRun),[1 2],'omitnan');
    meanCorrB(iRun) = median(withinB(chB(1):chB(2),chB(1):chB(2),iRun),[1 2],'omitnan');
end

withinAT = reshape(withinA,[31*31*4 1]); 
withinBT = reshape(withinB,[30*30*4 1]); 
% withinAB = reshape(mean([withinAT withinBT],2,'omitnan'),[31 31 4]); 


%%
titleLabel = {'1.2%','2.25%','3.25%', '1%' };
figure; 
for iFig = 1:4
    subplot(3,4,iFig); 
    imagesc(imgaussfilt(squeeze(btwAB(:,:,iFig)))); colormap parula; clim([0 0.7]); colorbar; 
    shading interp; axis square;
    title(['Btw probe: ' titleLabel{iFig}]); 

     subplot(3,4,iFig+4);
    imagesc(tril(imgaussfilt(squeeze(withinA(:,:,iFig))))); colormap parula; clim([0 0.7]); colorbar;
    shading interp; axis square; xticks(0:5:30); yticks(0:5:30);
    title(['Wthn probe: ' titleLabel{iFig}]);

    subplot(3,4,iFig+8);
    if iFig~=1
        x1 = smoothdata2(tril(squeeze(withinB(:,:,iFig)),-1));
        x1(logical(eye(length(x1)))) = 1;
    else
        x1 = imgaussfilt(squeeze(withinB(:,:,iFig)));
    end
    imagesc(tril(x1)); colormap parula; clim([0 0.7]); colorbar;
    shading interp; axis square;xticks(0:5:30); yticks(0:5:30);
    title(['Wthn probe: ' titleLabel{iFig}]);
end


%% Get power
maxRuns = max(cell2mat(cellfun(@length, datFileNumAll,'un',0)));
maxDates = size(allDates,1);

params.Fs     = fs;
params.fpass  = [6 120];
% params.pad    = 0;
params.tapers = [2 3];
psdA     = cellfun(@(x) NaN(119538,32,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);
psdB =  cellfun(@(x) NaN(119538,32,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);

tic;
for iDate = 1: size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        clear probeA probeB chA chB
        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);
        clear probeA probeB rawA rawB

        probeA = allProbeData{fileNum,iDate}.probe1Ch;
        probeB = allProbeData{fileNum,iDate}.probe2Ch;

        % rawA   = allProbeData{fileNum,iDate}.raw1Ch;
        % rawB   = allProbeData{fileNum,iDate}.raw2Ch;

        % Remove bad channels
        probeA(:,badElecA{fileNum,iDate}) = [];
        probeB(:,badElecB{fileNum,iDate}) = [];

        % rawA(:,badElecA{fileNum,iDate}) = [];
        % rawB(:,badElecB{fileNum,iDate}) = [];

        % Remove bad times
        probeA(allBadTimes{fileNum,iDate},:) = [];
        probeB(allBadTimes{fileNum,iDate},:) = [];
        %
        % rawA(allBadTimes{fileNum,iDate},:) = [];
        % rawB(allBadTimes{fileNum,iDate},:) = [];

        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);

        if chA(1) == 0 || chB(1) == 0 || isempty(probeA) || isempty(probeB)
            continue;
        end

        [psdA{iDate,fileNum},freq{iDate,fileNum}] = calculatePSD(probeA,params);
        psdB{iDate,fileNum} = calculatePSD(probeB,params);
    end
end

toc;

%% 
singleChRow = cellfun(@(x) isscalar(x),allVars.intraCorrBR(:,1));
freqLen = cell2mat(reshape(cellfun(@length,freq,'UniformOutput',0),[],1)); 
freqLen(allVars.removeDataIdx) = []; freqLen(singleChRow) = [];
fLen = unique(freqLen); fIdx = freqLen==fLen(3); 
avgPSDA = reshape(cellfun(@(x) median(x,2,'omitnan'),psdA,'UniformOutput',0),[],1);
avgPSDA(allVars.removeDataIdx) = []; avgPSDA(singleChRow) = []; 
avgPSDA(~fIdx) = []; avgPSDA = cat(2,avgPSDA{:});

avgPSDB = reshape(cellfun(@(x) median(x,2,'omitnan'),psdB,'UniformOutput',0),[],1);
avgPSDB(allVars.removeDataIdx) = []; avgPSDB(singleChRow) = []; 
avgPSDB(~fIdx) = []; avgPSDB = cat(2,avgPSDB{:});

allPSD = [avgPSDA avgPSDB];
medPSD = median(allPSD,2,'omitnan');
semVal = mad(allPSD,1,2)/sqrt(size(allPSD,2));

figure; subplot(121); plot(freq{1,1},10.*log10(median(avgPSDA,2))); set(gca, 'XScale', 'log');
ylim([-30 30]); ylabel('Power (dB)'); xlabel('Frequency (Hz)'); 
grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on'); box off; axis square;
subplot(122); plot(freq{1,1},10.*log10(median(avgPSDB,2))); set(gca, 'XScale', 'log');
ylim([-30 30]); ylabel('Power (dB)'); xlabel('Frequency (Hz)'); 
grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on'); box off; axis square; 
%%
l = 1:size(freq{1,1},2);
figure; 
plot(freq{1,1}(l),smooth(10.*log10(medPSD(l)),0.001)); hold on;
fill([freq{1,1}(l) fliplr(freq{1,1}(l))],...
    smooth(10.*log10([medPSD(l)+2.*semVal(l); flipud(medPSD(l)-2.*semVal(l))]),0.001),...
    'blue','FaceAlpha',0.3,'EdgeColor','none');
set(gca, 'XScale', 'log');
grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on'); box off; axis square; ylim([-25 15]);
 ylabel('Power (dB)'); xlabel('Frequency (Hz)'); 

 exportgraphics(gcf, 'C:\Users\kem294\OneDrive - University of Pittsburgh\Lab\Papers\Dual Probe Paper\Figures\MATLAB\Figure 6\avgPSD.eps', 'ContentType', 'vector', 'Resolution', 300);

 %%
 medPSD_smooth = smooth(10.*log10(medPSD), 0.001);
upper_smooth = smooth(10.*log10(medPSD + 2.*semVal), 0.001);
lower_smooth = smooth(10.*log10(medPSD - 2.*semVal), 0.001);

% Now downsample
downsample_factor = 100; % can be aggressive since already smoothed
l_down = l(1:downsample_factor:end);

figure;
plot(freq{1,1}(l_down), medPSD_smooth(l_down)); 
hold on;

fill([freq{1,1}(l_down) fliplr(freq{1,1}(l_down))],...
    [upper_smooth(l_down); flipud(lower_smooth(l_down))],...
    'blue','FaceAlpha',0.3,'EdgeColor','none');

set(gca, 'XScale', 'log');
grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on'); 
box off; axis square; ylim([-25 15]);
ylabel('Power (dB)'); xlabel('Frequency (Hz)');

%% Get coherence

params.Fs     = fs;
params.fpass  = [6 120];
% params.pad    = 0;
params.tapers = [2 3];
params.err    = [2 0.05];
params.trialave = 1;

paramsRaw       = params;
paramsRaw.fpass = [250 500];



winSize  = 1e3;
C     = cellfun(@(x) NaN(116,231,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);
confC = cellfun(@(x) NaN(1,231,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);
Cerr  = cellfun(@(x) NaN(2,116,231,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);

if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') 
    varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    varIdx  = ismember('allCohVals',varInfo);
else 
    varIdx = 1;
end

if varIdx
    allCohVals = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
        'allCohVals');
    allCohVals = allCohVals.allCohVals;
else
    allCohVals.coh = C;
    allCohVals.confThresh = confC;
    allCohVals.confIntJackknife = Cerr;

end

%% Get coherence in BLP
CBand     = cellfun(@(x) NaN(123,231,4,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);
confCBand = cellfun(@(x) NaN(1,231,4,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);
CerrBand  = cellfun(@(x) NaN(2,123,231,4,'single'), cell(maxDates,maxRuns), 'UniformOutput', false);
paramsNew = params;
paramsNew.fpass = [0 120]; 

for iDate = 1: size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        clear probeA probeB chA chB
        clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

        % if  all(isnan(allCohVals.coh{iDate,fileNum}(:,1)))%isempty(allCohVals.coh{iDate,iRun})
        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\coherenceVals_' num2str(fileNum) '.mat' ],'file') || ...
                all(isnan(allCohVals.coh{iDate,fileNum}(:,1))) 
            tic;
            clear probeA probeB rawA rawB
            probeA = allProbeData{fileNum,iDate}.probe1Ch;
            probeB = allProbeData{fileNum,iDate}.probe2Ch;

            % rawA   = allProbeData{fileNum,iDate}.raw1Ch;
            % rawB   = allProbeData{fileNum,iDate}.raw2Ch;

            % Remove bad channels
            probeA(:,badElecA{fileNum,iDate}) = [];
            probeB(:,badElecB{fileNum,iDate}) = [];

            % rawA(:,badElecA{fileNum,iDate}) = [];
            % rawB(:,badElecB{fileNum,iDate}) = [];

            % Remove bad times
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];
            %
            % rawA(allBadTimes{fileNum,iDate},:) = [];
            % rawB(allBadTimes{fileNum,iDate},:) = [];

            chA = estChInCortexA{iDate}(iRun,:);
            chB = estChInCortexB{iDate}(iRun,:);

            if chA(1) == 0 || chB(1) == 0 || isempty(probeA) || isempty(probeB)
                continue;
            end

            clear probeChACortex rawChACortex probeChBCortex rawChBCortex specB specA specBRaw specARaw
            probeChACortex = probeA(:,chA(1):chA(2));
            % rawChACortex   = rawA(:,chA(1):chA(2));

            probeChBCortex = probeB(:,chB(1):chB(2));
            % rawChBCortex   = rawA(:,chB(1):chB(2));

            dataLen = size(probeChBCortex,1);

            % Get all channel combinations
            combAllChT = table2array(combinations(1:size(probeChACortex,2),1:size(probeChBCortex,2)));
            if ~(chA(1)==1 || chB(2)==1)
                combAllCh{iDate,fileNum} = combAllChT(combAllChT(:,1)<=combAllChT(:,2),:);
            else
                combAllCh{iDate,fileNum} = combAllChT;
            end

            nWin     = 200;%numel(winStart);

            while winSize*nWin >dataLen % Check if the length of data exceeds the window size x number of windows
                nWin = nWin-1;
            end

            % Reshape data
            probeAR = reshape(probeChACortex(1:winSize*nWin,:), winSize, nWin, size(probeChACortex,2));
            probeBR = reshape(probeChBCortex(1:winSize*nWin,:), winSize, nWin, size(probeChBCortex,2));

            % Loop through channel pairs
            for iPair = 1:size(combAllCh{iDate,fileNum},1)
                [C{iDate,fileNum}(:,iPair),~,~,~,~,f,confC{iDate,fileNum}(:,iPair),~,Cerr{iDate,fileNum}(:,:,iPair)]=...
                    coherencyc(probeAR(:,:,combAllChT(iPair,1)),probeBR(:,:,combAllChT(iPair,2)),params);
            end           

            for iBand = 1:4
                switch iBand
                    case 1
                        xA = envelope(abs(filtfilt(bT,aT,double(probeChACortex))),5);
                        xB = envelope(abs(filtfilt(bT,aT,double(probeChBCortex))),5);
                    case 2
                        xA = envelope(abs(filtfilt(bA,aA,double(probeChACortex))),5);
                        xB = envelope(abs(filtfilt(bA,aA,double(probeChBCortex))),5);
                    case 3
                        xA = envelope(abs(filtfilt(bB,aB,double(probeChACortex))),5);
                        xB = envelope(abs(filtfilt(bB,aB,double(probeChBCortex))),5);
                    case 4
                        xA = envelope(abs(filtfilt(bG,aG,double(probeChACortex))),5);
                        xB = envelope(abs(filtfilt(bG,aG,double(probeChBCortex))),5);
                end

                xA = single(reshape(xA(1:winSize*nWin,:), winSize, nWin, size(probeChACortex,2)));
                xB = single(reshape(xB(1:winSize*nWin,:), winSize, nWin, size(probeChBCortex,2)));

                % Loop through channel pairs
                parfor iPair = 1:size(combAllCh{iDate,fileNum},1)
                    if iPair==1
                    [CBandT(:,iPair,iBand),~,~,~,~,fBandT(:,iPair),confCBandT(:,iPair,iBand),...
                        ~,CerrBandT(:,:,iPair,iBand)]=...
                        coherencyc(xA(:,:,combAllChT(iPair,1)),xB(:,:,combAllChT(iPair,2)),paramsNew);
                    else
                         [CBandT(:,iPair,iBand),~,~,~,~,~,confCBandT(:,iPair,iBand),...
                        ~,CerrBandT(:,:,iPair,iBand)]=...
                        coherencyc(xA(:,:,combAllChT(iPair,1)),xB(:,:,combAllChT(iPair,2)),paramsNew);
                    end
                end
            end

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\coherenceVals_' num2str(fileNum) '.mat' ],...
                'CBandT','confCBandT','CerrBandT','fBandT');

            CBand{iDate,fileNum}     = CBandT;
            confCBand{iDate,fileNum} = confCBandT;
            CerrBand{iDate,fileNum}  = CerrBandT;
            fBand{iDate,fileNum}     = fBandT;
            toc;

        else
            clear vals
            vals = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\coherenceVals_' num2str(fileNum) '.mat' ],...
                'CBandT','confCBandT','CerrBandT','fBandT');
            CBand{iDate,fileNum}     = vals.CBandT;
            confCBand{iDate,fileNum} = vals.confCBandT;
            CerrBand{iDate,fileNum}  = vals.CerrBandT;
            fBand{iDate,fileNum}     = vals.fBandT;

            C{iDate,fileNum}         = allCohVals.coh{iDate,fileNum};
            confC{iDate,fileNum}     = allCohVals.confThresh{iDate,fileNum};
            Cerr{iDate,fileNum}      = allCohVals.confIntJackknife{iDate,fileNum};
            combAllCh{iDate,fileNum} = allCohVals.combAllCh{iDate,fileNum} ;
            f                        = allCohVals.freqVals;
        end
    end
end

if ~varIdx
allCohVals.coh = C;
allCohVals.confThresh = confC;
allCohVals.confIntJackknife = Cerr;
save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
    'allCohVals','-append');
end

%% Visualize the mean coherence for BLP
fBandT = fBand{1,1}';
allCohBand = cellfun(@(x) squeeze(median(x,2,'omitnan')),CBand,'un',0);
allCohBand = cat(4,allCohBand{:});
allCohBand(:,:,allVars.removeDataIdx) = [];

allCerrBand = cellfun(@(x) squeeze(median(x,3,'omitnan')),CerrBand,'un',0);
allCerrBand = cat(4,allCerrBand{:}); allCerrBand(:,:,:,allVars.removeDataIdx) = [];
confCAvgBand = cellfun(@(x) squeeze(median(x,'all')),confCBand,'UniformOutput',0);
confCAvgBand = median(cat(2,confCAvgBand{:}),'all','omitnan');

save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
    'allCohBand','allCerrBand','confCAvgBand','fBandT','-append');

%%
figure;
for iBand = 1:4  
    subplot(2,2,iBand)
    plot(fBandT(2:end),median(squeeze(allCohBand((2:end),iBand,:)),2,'omitnan')); hold on;
    patch([fBandT(2:end) fliplr(fBandT(2:end))], [median(squeeze(allCerrBand(1,(2:end),iBand,:)),2,'omitnan');...
        flipud(median(squeeze(allCerrBand(2,(2:end),iBand,:)),2,'omitnan'))],'blue','FaceAlpha',0.3,'EdgeColor','none')
   set(gca, 'XScale', 'log'); 
    box off; xlabel('Frequency (Hz)'); ylabel('Coherence'); title(bandLabels{iBand});
    yline(confCAvgBand,'LineWidth',1); ylim([0 1]);
    grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on');

end

%% Plot coherence versus frequencies
connValsR = connValsAll';
connValsR = reshape(connValsR,[maxRuns*maxDates 1]);
connValsR(allVars.removeDataIdx) = [];

distValsR = distSitesAll';
distValsR = reshape(distValsR,[maxRuns*maxDates 1]);
distValsR(allVars.removeDataIdx) = [];

%%
cohNew  = cat(3,allCohBand,cohBandCS);
cErrNew = cat(4,allCerrBand,CerrBandCS);
connAll = [connValsR; connValsCS];
distAll = [distValsR; distValsCS];

%%
figure;
for iBand = 1:4  
    subplot(2,2,iBand)
    plot(fBandT(2:end),median(squeeze(cohNew((2:end),iBand,:)),2,'omitnan')); hold on;
    patch([fBandT(2:end) fliplr(fBandT(2:end))], [median(squeeze(cErrNew(1,(2:end),iBand,:)),2,'omitnan');...
        flipud(median(squeeze(cErrNew(2,(2:end),iBand,:)),2,'omitnan'))],'blue','FaceAlpha',0.3,'EdgeColor','none')
   set(gca, 'XScale', 'log'); 
    box off; xlabel('Frequency (Hz)'); ylabel('Coherence'); title(bandLabels{iBand});
    yline(confCAvgBand,'LineWidth',1); ylim([0 1]);
    grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on'); xlim([1 10]);xticks(1:10);

end

%% Separate the coherences based on bands 
clear cohValsFreqTrace
for iTrace = 1:4
    for iBand = 1:2
        switch iBand
            case 1
                freqIdx = fBandT>=1 & fBandT<=4;
            case 2
                freqIdx = fBandT>4 & fBandT<=8;
            % case 3
            %     freqIdx = fBandT>=betaBand(1) & fBandT<=betaBand(2);
            % case 4
            %     freqIdx = fBandT>=gammaBand(1) & fBandT<=gammaBand(2);
        end
        cohValsFreqTrace(:,iBand,iTrace) = median(squeeze(cohNew(freqIdx,iTrace,:)),1,'omitnan')';
    end
end
% cohValsFreqR  = reshape(cohValsFreq,[maxRuns*maxDates 4]);

% threshCond = sum(cohValsFreq<confCAvg,2)>=1;
% cohValsFreq(threshCond,:) = []; connAll(threshCond) = []; distAll(threshCond) = [];

%%
freqLabels = {'1-4 Hz'; '4-8 Hz'};
for iTrace = 1:4
    figure;
    for iBand = 1:2
        subplot(1,2,iBand);
        plot(connAll,cohValsFreqTrace(:,iBand,iTrace),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
        ylim([0 1]); axis square;
        xVal = connAll; yVal = cohValsFreqTrace(:,iBand,iTrace);
        coeff = fit(xVal,double(yVal),'poly1','Robust','LAR');
        xFit  = linspace(min(xVal),max(xVal),1000);
        yFit  = coeff.p1*xFit + coeff.p2; mdl = fitlm(xVal,yVal,'RobustOpts','on');
        plot(xFit,yFit,'-k','LineWidth',1);
        text(0.8, 0.35,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
        text(0.8, 0.3,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]);
        title(freqLabels{iBand}); xlabel('Functional connectivity');
        ylabel('Coherence');
    end
    sgtitle(['BLP: ' bandLabels{iTrace}]);
end

% y = zscore(connValsR);
% X  = zscore(cohValsFreqR);
% clear  beta sigma eVal covMat

% Get the linear model
% mdlAll= fitlm(cohValsFreq,connAll,'VarNames',{'Theta','Alpha','Beta','Gamma','FC'});

%%
for iTrace = 1:4
    figure;

for iBand = 1:2
    subplot(1,2,iBand);
    plot(distAll,cohValsFreqTrace(:,iBand,iTrace),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
    ylim([0 1]); axis square;
    xVal = distAll; yVal = cohValsFreqTrace(:,iBand,iTrace);

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
    
    text(8, 0.75,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(8, 0.7,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]); 
    title(freqLabels{iBand}); xlabel('Distance (mm)');
    ylabel('Coherence');
end
sgtitle(['BLP: ' bandLabels{iTrace}]);
end

%%
confC = cellfun(@(x) single(x),confC,'uniformoutput',0);
allCohVals.coh = C;
allCohVals.confThresh = confC;
allCohVals.confIntJackknife = Cerr;
allCohVals.freqVals = f;
% allCohVals.combAllCh = combAllCh;
% save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
%     'allCohVals','-append');


%% Visualize mean coherence of all recordings
figure; idx =1;
for iDate = 1:size(C,1)
    for iRun = 1:size(C,2)
        if isempty(C{iDate,iRun})| isnan(C{iDate,iRun}); idx = idx+1; continue; end
        subplot(size(C,1), size(C,2), idx); plot(f,median(C{iDate,iRun},2,'omitnan')); hold on;
       patch([f fliplr(f)], [median(squeeze(Cerr{iDate,iRun}(1,:,:)),2,'omitnan');...
            flipud(median(squeeze(Cerr{iDate,iRun}(2,:,:)),2,'omitnan'))],'blue','FaceAlpha',0.3,'EdgeColor','none')
        box off; xlabel(' Frequency (Hz)'); ylabel('Coherence');
        yline(median(confC{iDate,iRun},"all"),'LineWidth',1); ylim([0 0.5]);
        idx = idx+1; set(gca, 'XScale', 'log');
    end
end

%% Plot the median coherence
allCohRec = cell2mat(reshape(cellfun(@(x) (median(x,2,'omitnan')),C,'un',0),[size(C,1)*size(C,2),1])');
allCohRec(:,allVars.removeDataIdx) = [];
allCerr = cellfun(@(x) (median(x,3,'omitnan')),Cerr,'un',0);
allCerr = cat(3,allCerr{:}); allCerr(:,:,allVars.removeDataIdx) = [];
confCAvg = median(cell2mat(cellfun(@(x) median(x,2),confC,'UniformOutput',0)),'all','omitnan');

l = load(['D:\Data\CharlieSheen_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
    'allCohRec','allCerr','confCAvg');

cohCS = l.allCohRec;
allCerrCS = l.allCerr;

% save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
%     'allCohRec','allCerr','confCAvg','-append');
%%
figure;
plot(f,median(allCohRec,2,'omitnan')); hold on;
  patch([f fliplr(f)], [median(squeeze(allCerr(1,:,:)),2,'omitnan');...
            flipud(median(squeeze(allCerr(2,:,:)),2,'omitnan'))],'blue','FaceAlpha',0.3,'EdgeColor','none')
  box off; xlabel(' Frequency (Hz)'); ylabel('Coherence');
  yline(confCAvg,'LineWidth',1); ylim([0 0.5]); set(gca, 'XScale', 'log');
 grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on');


 %%
 cohNew = [allCohRec cohCS];
 cErrNew = cat(3,allCerr,allCerrCS);
 connAll = [connValsR; connValsCS];
 distAll = [distValsR; distValsCS];
%%
figure;
plot(f,median(cohNew,2,'omitnan')); hold on;
  patch([f fliplr(f)], [median(squeeze(cErrNew(1,:,:)),2,'omitnan');...
            flipud(median(squeeze(cErrNew(2,:,:)),2,'omitnan'))],'blue','FaceAlpha',0.3,'EdgeColor','none')
  box off; xlabel(' Frequency (Hz)'); ylabel('Coherence');
  yline(confCAvg,'LineWidth',1); ylim([0 0.5]); set(gca, 'XScale', 'log');
 grid on; set(gca, 'YGrid', 'off', 'XGrid', 'on');



%% Separate the coherences based on bands 
for iBand = 1:4
    switch iBand
        case 1
            freqIdx = f>=thetaBand(1) & f<=thetaBand(2);
        case 2
            freqIdx = f>=alphaBand(1) & f<=alphaBand(2);
        case 3
            freqIdx = f>=betaBand(1) & f<=betaBand(2);
        case 4
            freqIdx = f>=gammaBand(1) & f<=gammaBand(2);
    end
    cohValsFreq(:,iBand) = median(cohNew(freqIdx,:),1,'omitnan')';
end

% cohValsFreqR  = reshape(cohValsFreq,[maxRuns*maxDates 4]);

threshCond = sum(cohValsFreq<confCAvg,2)>=1;
cohValsFreq(threshCond,:) = []; connAll(threshCond) = []; distAll(threshCond) = [];

%%
figure;
for iBand = 1:4
    subplot(2,2,iBand);
    plot(connAll,cohValsFreq(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
    ylim([0 0.8]); axis square;
    xVal = connAll; yVal = cohValsFreq(:,iBand);
    coeff = fit(xVal,double(yVal),'poly1','Robust','LAR');
    xFit  = linspace(min(xVal),max(xVal),1000);
    yFit  = coeff.p1*xFit + coeff.p2; mdl = fitlm(xVal,yVal,'RobustOpts','on');
    plot(xFit,yFit,'-k','LineWidth',1);
    text(0.8, 0.35,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(0.8, 0.3,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]); 
    title(bandLabels{iBand}); xlabel('Functional connectivity');
    ylabel('Coherence');
end

% y = zscore(connValsR);
% X  = zscore(cohValsFreqR);
% clear  beta sigma eVal covMat

% Get the linear model
mdlAll= fitlm(cohValsFreq,connAll,'VarNames',{'Theta','Alpha','Beta','Gamma','FC'});

%%
figure;
for iBand = 1:4
    subplot(2,2,iBand);
    plot(distAll,cohValsFreq(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
    ylim([0 0.8]); axis square;
    xVal = distAll; yVal = cohValsFreq(:,iBand);

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

    
    text(8, 0.75,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(8, 0.7,['p-val: ' num2str(mdl.ModelFitVsNullModel.Pvalue)]); 
    title(bandLabels{iBand}); xlabel('Distance (mm)');
    ylabel('Coherence');
end

%% Divide coherence into three compartments
chSplit = 6;
superCh = 1:chSplit;
midCh   = chSplit+1 : 2*chSplit;
cohDepth = NaN(maxDates,maxRuns,9,116,'single');

for iDate = 1:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    
    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);
        
        if chA(1) == 0 || chB(1) == 0; continue; end
        
        deepChA = 2*chSplit+1:max(combAllCh{iDate,fileNum}(:,1));
        deepChB = 2*chSplit+1:max(combAllCh{iDate,fileNum}(:,2));
        
        % Original combinations
        comb = combAllCh{iDate,fileNum};
        C_vals = C{iDate,fileNum};
        
        % Create reverse pairs to handle symmetry
        comb_reverse = [comb(:,2), comb(:,1)];
        comb_full = [comb; comb_reverse];
        C_full = [C_vals, C_vals];  % Duplicate coherence values (along frequency dimension)
        
        % Now calculate cohDepth with bidirectional mask
        % Format: mask1 | mask2 (check both directions)
        
        % 1. Super A ↔ Super B
        mask = (ismember(comb_full(:,1), superCh) & ismember(comb_full(:,2), superCh));
        cohDepth(iDate,fileNum,1,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 2. Super A ↔ Mid B
        mask = (ismember(comb_full(:,1), superCh) & ismember(comb_full(:,2), midCh));
        cohDepth(iDate,fileNum,2,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 3. Super A ↔ Deep B
        mask = (ismember(comb_full(:,1), superCh) & ismember(comb_full(:,2), deepChB));
        cohDepth(iDate,fileNum,3,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 4. Mid A ↔ Super B 
        mask = (ismember(comb_full(:,1), midCh) & ismember(comb_full(:,2), superCh));
        cohDepth(iDate,fileNum,4,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 5. Mid A ↔ Mid B
        mask = (ismember(comb_full(:,1), midCh) & ismember(comb_full(:,2), midCh));
        cohDepth(iDate,fileNum,5,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 6. Mid A ↔ Deep B
        mask = (ismember(comb_full(:,1), midCh) & ismember(comb_full(:,2), deepChB));
        cohDepth(iDate,fileNum,6,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 7. Deep A ↔ Super B (now will find pairs!)
        mask = (ismember(comb_full(:,1), deepChA) & ismember(comb_full(:,2), superCh));
        cohDepth(iDate,fileNum,7,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 8. Deep A ↔ Mid B (now will find pairs!)
        mask = (ismember(comb_full(:,1), deepChA) & ismember(comb_full(:,2), midCh));
        cohDepth(iDate,fileNum,8,:) = median(C_full(:, mask), 2, 'omitnan');
        
        % 9. Deep A ↔ Deep B
        mask = (ismember(comb_full(:,1), deepChA) & ismember(comb_full(:,2), deepChB));
        cohDepth(iDate,fileNum,9,:) = median(C_full(:, mask), 2, 'omitnan');
    end
end

cohDepthR = reshape(cohDepth, [maxRuns*maxDates 3 3 116]);
cohDepthR(allVars.removeDataIdx,:,:,:) = [];
% cohDepthR(nanValsCoh,:,:,:) = [];
cohDepthR(threshCond,:,:,:) = [];

%% Plot the heatmap for all frequencies
figure;
for iBand = 1:4
    switch iBand
        case 1
            freqIdx = f>=thetaBand(1) & f<=thetaBand(2);
        case 2
            freqIdx = f>=alphaBand(1) & f<=alphaBand(2);
        case 3
            freqIdx = f>=betaBand(1) & f<=betaBand(2);
        case 4
            freqIdx = f>=gammaBand(1) & f<=gammaBand(2);
    end
    cohValsFreqDepth(:,:,iBand) = squeeze(median(cohDepthR(:,:,:,freqIdx),[1 4],'omitnan'));
    subplot(2,2,iBand); imagesc(cohValsFreqDepth(:,:,iBand)); 
    if iBand<=3,clim([0 0.45]); else, clim([0 0.2]); end
    title(bandLabels{iBand});
    xlabel('Probe B'); ylabel('Probe A'); axis square; colorbar;
    xticks(1:3); yticks(1:3); xticklabels({'S','M','D'});yticklabels({'S','M','D'});
end


%% Perform Phase amplitude coupling
% Get the PAC comodulogram
% Set the amplitude range
clear allPACVars 
gammaRange = gammaBand(1):5:gammaBand(2);
gammaRange(gammaRange>=60 & gammaRange<=65)=[];

% Get the phase range
lowFreqRange = 6:2:30;

nHigh    = size(gammaRange,2);
nLow     = size(lowFreqRange,2);

% Initialize window and step size
winSize  = 10e3;
stepSize = 10e3;

chSplit = 6;

% Get the phase amplitude coallMupling for the recordings...
for iDate = 1:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    
    if strcmp(monkeyName,'Whiskey') && iDate == 1; continue; end
    for iRun = 1:length(datFileNum)

        fileNum = datFileNum(iRun);
    clear modIdxA2A modIdxB2B modIdxAamp2Bphase modIdxAamp2Bphase
        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);

        if chA(1)== 0 || chB(1)==0; continue; end

        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_10sec_' num2str(fileNum) '.mat'],'file') %
            clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);
            tic;

            [amplitudeA, amplitudeB, phaseA, phaseB] = calculatePhaseAmpSignals(monkeyName,expDate,hemisphere,fileNum,...
                allProbeData{fileNum,iDate}.probe1Ch,allProbeData{fileNum,iDate}.probe2Ch,badElecA{fileNum,iDate},...
                badElecB{fileNum,iDate},allBadTimes{fileNum,iDate},estChInCortexA{iDate}(iRun,:),...
                estChInCortexB{iDate}(iRun,:),gammaRange,lowFreqRange);

            % Get all channel combinations
            combAllCh    = table2array(combinations(1:size(amplitudeA,3),1:size(amplitudeB,3)));
            combAllChA2A = table2array(combinations(1:size(amplitudeA,3),1:size(amplitudeA,3)));
            combAllChB2B = table2array(combinations(1:size(amplitudeB,3),1:size(amplitudeB,3)));
           
            dataLen   = size(amplitudeA,2);

            % Calculate windows
            winStart = 1:stepSize:(dataLen-winSize);
            winEnd   = winStart+winSize-1;
            nWin     = numel(winStart);
            if nWin>60; nWin = 60; winStart = winStart(1:nWin); winEnd = winEnd(1:nWin); end

            %         [modIdxA2A,modIdxB2B,modIdxAamp2Bphase,modIdxAphase2Bamp] = getPhaseAmpCoupling(phaseA,...
            % amplitudeA,phaseB,amplitudeB,winStart,winEnd,winSize,combAllCh(:,1), combAllCh(:,2));

            modIdxA2A = getPhaseAmpCoupling_v2(phaseA,amplitudeA,winStart,winEnd,winSize,combAllChA2A(:,1),combAllChA2A(:,2));
            modIdxB2B = getPhaseAmpCoupling_v2(phaseB,amplitudeB,winStart,winEnd,winSize,combAllChB2B(:,1),combAllChB2B(:,2));

            modIdxAamp2Bphase = getPhaseAmpCoupling_v2(phaseB,amplitudeA,winStart,winEnd,winSize,combAllCh(:,1), combAllCh(:,2));
            modIdxAphase2Bamp = getPhaseAmpCoupling_v2(phaseA,amplitudeB,winStart,winEnd,winSize,combAllCh(:,2),combAllCh(:,1)); 

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_10sec_' num2str(fileNum) '.mat'],...
                'modIdxA2A','modIdxB2B','modIdxAamp2Bphase','modIdxAphase2Bamp','combAllCh','combAllChA2A','combAllChB2B');
        toc;
        end

        allPACVars{fileNum,iDate} = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramVals_10sec_' num2str(fileNum) '.mat']);
    end
end


%% Getting the surrogate/shuffled distribution
clear allCtrlVals
nHigh = size(gammaRange,2);
nLow  = size(lowFreqRange,2);

winSize  = 10e3;
stepSize = 10e3;

% shiftLen = [1 5 10 20 50 100].*1e3;
% nShift   = length(shiftLen);

for iDate = 1:size(allDates,1)
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    if strcmp(monkeyName,'Whiskey') && iDate == 1; continue; end
    for iRun = 1:length(datFileNum)
        % if iDate == 3 && iRun>8; continue; end 

        fileNum = datFileNum(iRun);
        chA = estChInCortexA{iDate}(iRun,:);
        chB = estChInCortexB{iDate}(iRun,:);
        if chA(1)== 0 || chB(1)==0; continue; end
        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate...
                '\Electrophysiology\modulogramCtrl10sec_' num2str(fileNum) '.mat'],'file')
            clc; disp(['Getting shuffled/surrogate data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);
  

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

            winStart = 1:stepSize:(dataLen-winSize);
            winEnd   = winStart+stepSize-1;
            nWin     = numel(winStart);
            if nWin>60; nWin = 60; winStart = winStart(1:nWin); winEnd = winEnd(1:nWin); end

            % Create artificial surrogates for the data
            % Shuffle every sample of the phase time course
            clear phaseAShuff phaseBShuff
            phaseAShuff = zeros([size(phaseA)],'single');
            phaseBShuff = zeros([size(phaseB)],'single');
           
            rng('shuffle');
            comb1 = randperm(dataLen);
            phaseAShuff = phaseA(:,comb1,:);
            phaseBShuff = phaseB(:,comb1,:); 

            tic;

            modIdxAllA2AShuffleT = getPhaseAmpCoupling_v2(phaseAShuff,amplitudeA,winStart,winEnd,winSize,(1:nChan)',(1:nChan)');
            modIdxAllB2BShuffleT = getPhaseAmpCoupling_v2(phaseBShuff,amplitudeB,winStart,winEnd,winSize,(1:nChan)',(1:nChan)');

            modIdxAllA2BShuffleT = getPhaseAmpCoupling_v2(phaseBShuff,amplitudeA,winStart,winEnd,winSize,(1:nChan)',(1:nChan)');
            modIdxAllB2AShuffleT = getPhaseAmpCoupling_v2(phaseAShuff,amplitudeB,winStart,winEnd,winSize,(1:nChan)',(1:nChan)');
            toc;

            save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramCtrl10sec_' num2str(fileNum) '.mat'],...
                'modIdxAllA2BShuffleT','modIdxAllB2AShuffleT','modIdxAllA2AShuffleT','modIdxAllB2BShuffleT');
        end

            allCtrlVals{fileNum,iDate} = matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\modulogramCtrl10sec_' num2str(fileNum) '.mat']);
    end
end

% Organizing the control/surrogate distribution
zeroVals = reshape(cell2mat(cellfun(@(x) isempty(x),allCtrlVals,'un',0))',[numel(allCtrlVals) 1]) | allVars.removeDataIdx;
allCtrlValsT = reshape(allCtrlVals',[numel(allCtrlVals) 1]);
allCtrlValsT(zeroVals) = [];

for iL = 1:size(allCtrlValsT,1)
    a2bShuffle(iL,:) = reshape(squeeze(median(allCtrlValsT{iL}.modIdxAllA2BShuffleT,[1 4],'omitnan')),[nLow*nHigh,1]);
    b2aShuffle(iL,:) = reshape(squeeze(median(allCtrlValsT{iL}.modIdxAllB2AShuffleT,[1,4],'omitnan')),[nLow*nHigh,1]); 
    a2aShuffle(iL,:) = reshape(squeeze(median(allCtrlValsT{iL}.modIdxAllA2AShuffleT,[1,4],'omitnan')),[nLow*nHigh,1]);
    b2bShuffle(iL,:) = reshape(squeeze(median(allCtrlValsT{iL}.modIdxAllB2BShuffleT,[1 4],'omitnan')),[nLow*nHigh,1]); 
end

% Getting the 99th percentile value and setting it as threshold
a2bThresh = prctile(a2bShuffle,99)';
b2aThresh = prctile(b2aShuffle,99)';
a2aThresh = prctile(a2aShuffle,99)';
b2bThresh = prctile(b2bShuffle,99)'; 

%% Plotting an example
figure; 
for iDate = 4%2:7
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    for iRun = 3%1:length(datFileNum)
        if isempty(allPACVars{iRun,iDate}); continue; end
         tempAvgBtwMI = reshape(squeeze(median(cat(3,squeeze(median(allPACVars{iRun,iDate}.modIdxAamp2Bphase,[1 4],'omitnan')),...
            squeeze(median(allPACVars{iRun,iDate}.modIdxAphase2Bamp,[1 4],'omitnan'))),3,'omitnan')),[nHigh*nLow 1]);
         tempAvgBtwMI(tempAvgBtwMI<a2bThresh | tempAvgBtwMI <b2aThresh) = 0;
         % contourf(reshape(tempAvgBtwMI,[nHigh nLow]),'lines','none'); 
         imagesc(imgaussian(reshape(tempAvgBtwMI,[nHigh nLow]),1));
         set(gca,'YDir','normal'); clim([0 5e-4]);colorbar; axis square;shading interp;
         % set(gca,'YDir','normal'); clim([0 6e-4]);colorbar; axis square;
         title(['date: ' expDate ' run:' num2str(iRun)]); yticks(1:nHigh); xticks(1:nLow);
         yticklabels(gammaRange); xticklabels(lowFreqRange);
    end
end

%% Synthesizing the results
 clear thetaHighGammaPairsBtw thetaHighGammaPairsWithin thetaLowGammaPairsWithin thetaLowGammaPairsBtw...
     avgMIWithin avgMIBetween pacVals

minVal =  -1;
maxVal = 1;

[phaseVal, ampVal] = meshgrid(lowFreqRange,gammaRange);
allCombVec = [phaseVal(:) ampVal(:)];

thetaHighGammaIdx = ismember(allCombVec,...
    single(table2array(combinations(lowFreqRange(lowFreqRange<=8),gammaRange(gammaRange>=70)))),'rows');

thetaLowGammaIdx = ismember(allCombVec,...
    single(table2array(combinations(lowFreqRange(lowFreqRange<8),gammaRange(gammaRange>30 & gammaRange<=50)))),'rows');

ctrlIdx    = ismember(allCombVec,...
    single(table2array(combinations(lowFreqRange(lowFreqRange>24 & lowFreqRange<30),gammaRange(gammaRange>50 & gammaRange<80)))),'rows');

for iDate = 1:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        if isempty(allPACVars{fileNum,iDate}); continue; end 

        % Between probes
        tempAvgBtwMI = reshape(squeeze(median(cat(3,squeeze(median(allPACVars{fileNum,iDate}.modIdxAamp2Bphase,[1 4],'omitnan')),...
            squeeze(median(allPACVars{fileNum,iDate}.modIdxAphase2Bamp,[1 4],'omitnan'))),3,'omitnan')),[nHigh*nLow 1]);
        
        % Within probes
        tempAvgWithinMI = reshape(squeeze(median(cat(3,squeeze(median(allPACVars{fileNum,iDate}.modIdxA2A,[1 4],'omitnan')),...
            squeeze(median(allPACVars{fileNum,iDate}.modIdxB2B,[1 4],'omitnan'))),3,'omitnan')),[nHigh*nLow 1]);

        % rankAvgBtw = tiedrank(tempAvgBtwMI); rankAvgBtw = minVal + ((rankAvgBtw-1)*(maxVal-minVal)/(length(rankAvgBtw)-1));
        % rankAvgWth = tiedrank(tempAvgWithinMI); rankAvgWth = minVal+ ((rankAvgWth-1)*(maxVal-minVal)/length(rankAvgWth)-1);

        ctrlMIBtw(fileNum,iDate)    = median(tempAvgBtwMI(ctrlIdx),'omitnan');
        ctrlMIWithin(fileNum,iDate) = median(tempAvgWithinMI(ctrlIdx),'omitnan');

        % Removing the pairs that do not have significant modulation
        tempAvgBtwMI(tempAvgBtwMI<a2bThresh | tempAvgBtwMI <b2aThresh) = NaN;
        tempAvgWithinMI(tempAvgWithinMI<a2aThresh | tempAvgWithinMI<b2bThresh) = NaN;

        avgMIBetween(fileNum,iDate,:) = tempAvgBtwMI;
        avgMIWithin(fileNum,iDate,:)  = tempAvgWithinMI;

        thetaLowGammaPairsBtw(fileNum,iDate)  = median(avgMIBetween(fileNum,iDate,thetaLowGammaIdx),'omitnan');
        thetaHighGammaPairsBtw(fileNum,iDate) = median(avgMIBetween(fileNum,iDate,thetaHighGammaIdx),'omitnan');

        thetaLowGammaPairsWithin(fileNum,iDate)  = median(avgMIWithin(fileNum,iDate,thetaLowGammaIdx),'omitnan');
        thetaHighGammaPairsWithin(fileNum,iDate) = median(avgMIWithin(fileNum,iDate,thetaHighGammaIdx),'omitnan');
        
    end

end


%%  Rearranging the variables
miMatSize     = numel(thetaHighGammaPairsWithin);
matSize       = size(avgMIBetween);
avgMIBetweenT = reshape(permute(avgMIBetween,[2 1 3]), [miMatSize nHigh*nLow]);
avgMIWithinT  = reshape(permute(avgMIWithin,[2 1 3]),[miMatSize nHigh*nLow]);

% zeroIdx                  = (all(avgMIBetweenT==0,2));
avgMIBetweenT(zeroVals,:) = [];
avgMIWithinT(zeroVals,:)  = [];

thetaHighGammaPairsBtw(thetaHighGammaPairsBtw==0) = NaN; 
thetaLowGammaPairsBtw(thetaLowGammaPairsBtw ==0)  = NaN; 

thetaHighGammaPairsWithin(thetaHighGammaPairsWithin==0) = NaN; 
thetaLowGammaPairsWithin(thetaLowGammaPairsWithin ==0)  = NaN; 

ctrlMIBtw(ctrlMIBtw==0) = NaN; 
ctrlMIWithin(ctrlMIWithin==0) = NaN;

thetaHGBtwT = reshape(thetaHighGammaPairsBtw',[miMatSize 1]);    thetaHGBtwT(zeroVals) = [];
thetaLGBtwT = reshape(thetaLowGammaPairsBtw',[miMatSize 1]);     thetaLGBtwT(zeroVals) = [];
thetaHGWthT = reshape(thetaHighGammaPairsWithin',[miMatSize 1]); thetaHGWthT(zeroVals) = [];
thetaLGWthT = reshape(thetaLowGammaPairsWithin',[miMatSize 1]);  thetaLGWthT(zeroVals) = [];

ctrlMIBtwT = reshape(ctrlMIBtw',[miMatSize 1]); ctrlMIBtwT(zeroVals)=[]; 
ctrlMIWithinT = reshape(ctrlMIWithin',[miMatSize 1]); ctrlMIWithinT(zeroVals) = [];


if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\allPACVars.mat'],'file')
    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\allPACVars.mat'],...
        'allPACVars','allCtrlVals','avgMIBetweenT','avgMIWithinT','thetaHGBtwT',...
        'thetaLGBtwT','thetaHGWthT','thetaLGWthT',...
        'ctrlMIBtwT','ctrlMIWithinT','a2bThresh','b2aThresh','a2aThresh','b2bThresh');

else
        save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\allPACVars.mat'],...
        'allPACVars','allCtrlVals','avgMIBetweenT','avgMIWithinT','thetaHGBtwT',...
        'thetaLGBtwT','thetaHGWthT','thetaLGWthT',...
        'ctrlMIBtwT','ctrlMIWithinT','a2bThresh','b2aThresh','a2aThresh','b2bThresh','-append');

end

%% Plotting the average 
figure;
miBtw = avgMIBetweenT;miBtw(isnan(miBtw))=0;
imagesc(imgaussian(reshape(median(miBtw,1),[nHigh nLow]),1));
% contourf(reshape(median(miBtw,1),[nHigh nLow]),'lines','none');
set(gca,'YDir','normal'); clim([0 5e-4]);colorbar; axis square;shading interp
% set(gca,'YDir','normal'); clim([0 6e-4]);colorbar; axis square;
yticks(1:nHigh); xticks(1:nLow); yticklabels(gammaRange); xticklabels(lowFreqRange);


%% Showing the distribution of modulation indices
figure; 
subplot(121); boxplot([median(avgMIBetweenT,2,'omitnan') median(avgMIWithinT,2,'omitnan')],{'Between probes','Within probes'}); 
box off; axis square; ylim([2.7e-4 4.3e-4]); ylabel('Average modulation index');
subplot(122);violin([median(avgMIBetweenT,2,'omitnan') median(avgMIWithinT,2,'omitnan')],{'Between probes','Within probes'}); 
box off; axis square; ylim([2.5e-4 4.5e-4]); ylabel('Average modulation index');

%% 
figure;
subplot(121); boxplot([thetaLGBtwT thetaHGBtwT ctrlMIBtwT thetaLGWthT thetaHGWthT ctrlMIWithinT],...
    {'Theta_LowGamma Between', 'Theta_highGamma Between', 'Ctrl Between','Theta_LowGamma Within',...
    'Theta_highGamma Within','Ctrl Within'});

box off; axis square; ylim([2.4e-4 7e-4]);ylabel('Average modulation index');

subplot(122); violin([thetaLGBtwT thetaHGBtwT ctrlMIBtwT thetaLGWthT thetaHGWthT ctrlMIWithinT],...
    {'Theta_LowGamma Between', 'Theta_highGamma Between', 'Ctrl Between','Theta_LowGamma Within',...
    'Theta_highGamma Within','Ctrl Within'});
box off; axis square; ylim([2.4e-4 8e-4]);ylabel('Average modulation index');

%% 
singleChRow = cellfun(@(x) isscalar(x),allVars.intraCorrBR(:,1));

connValsNew = allVars.connValsR; connValsNew(singleChRow) = [];
% connValsAll; connValsNew = reshape(connValsNew,[miMatSize 1]);
% connValsNew(zeroIdx) = [];

distValsNew = allVars.distValsR; distValsNew(singleChRow) = [];
% distSitesAll; distValsNew = reshape(distValsNew,[miMatSize 1]); 
% distValsNew(zeroIdx)= [];
%%
medMIBtw = median(avgMIBetweenT,2,'omitnan');
idxMI = tiedrank(medMIBtw);
idxMI = minVal + ((idxMI-1)*(maxVal-minVal)/(length(idxMI)-1));

lowGammaRank = tiedrank(thetaLGBtwT);
lowGammaRank = minVal + ((lowGammaRank-1)*(maxVal-minVal)/(length(lowGammaRank)-1));
highGammaRank = tiedrank(thetaHGBtwT);
highGammaRank = minVal + ((highGammaRank-1)*(maxVal-minVal)/(length(highGammaRank)-1));

ctrlRank      = tiedrank(ctrlMIBtwT);
ctrlRank      = minVal + ((ctrlRank-1)*(maxVal-minVal)/(length(ctrlRank)-1));
%%
figure;subplot(131);showLinearFit(connValsNew,lowGammaRank); axis square; title('Theta-low gamma');
subplot(132); showLinearFit(connValsNew,highGammaRank); axis square; title('Theta-high gamma');
subplot(133); showLinearFit(connValsNew,ctrlRank); axis square; title('Control');

%%
figure;subplot(121); showLinearFit(connValsNew,idxMI);axis square;
subplot(122);showLinearFit(connValsNew,medMIBtw);axis square;


%%
figure; 
showLinearFit(connValsNew,median(avgMIBetweenT,2,'omitnan')); 
ylim([2.8e-4 3.8e-4]); axis square;

%%
figure; 
subplot(121); 
showLinearFit(connValsNew,thetaLGBtwT); 
ylim([2.8e-4 7e-4]); axis square;
ylabel('Theta-low gamma MI');

subplot(122); 
showLinearFit(connValsNew,thetaHGBtwT); 
ylim([2.8e-4 7e-4]); axis square; 
xlabel('Functional connectivity'); 
ylabel('Theta-high gamma MI');

%% Laminar PAC

[phaseVal, ampVal] = meshgrid(lowFreqRange,gammaRange);
allCombVec = [phaseVal(:) ampVal(:)];

thetaHighGammaIdx = ismember(allCombVec,...
    single(table2array(combinations(lowFreqRange(lowFreqRange<=8),gammaRange(gammaRange>=70)))),'rows');

thetaLowGammaIdx = ismember(allCombVec,...
    single(table2array(combinations(lowFreqRange(lowFreqRange<8),gammaRange(gammaRange>30 & gammaRange<=50)))),'rows');

for iDate = 1:size(allDates,1)
    clear expDate datFileNum 
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};

    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);
        if isempty(allPACVars{fileNum,iDate}); continue; end 
               
        pacVals(fileNum,iDate).aAmpBPhase = allPACVars{fileNum,iDate}.modIdxAamp2Bphase;
        pacVals(fileNum,iDate).aPhaseBAmp = allPACVars{fileNum,iDate}.modIdxAphase2Bamp;
        pacVals(fileNum,iDate).a2a        = allPACVars{fileNum,iDate}.modIdxA2A; 
        pacVals(fileNum,iDate).b2b        = allPACVars{fileNum,iDate}.modIdxB2B; 

        chComb = allPACVars{fileNum,iDate}.combAllCh; 
        
        maxChA = max(chComb(:,1));
        maxChB = max(chComb(:,2));   

        aSuperbSuper = ismember(chComb,table2array(combinations(1:chSplit,1:chSplit)),'rows');
        aSuperbMid   = ismember(chComb,table2array(combinations(1:chSplit,chSplit+1:2*chSplit)),'rows');
        aSuperbDeep  = ismember(chComb,table2array(combinations(1:chSplit,2*chSplit+1:maxChB)),'rows');

        aMidbSuper = ismember(chComb,table2array(combinations(chSplit+1:2*chSplit,1:chSplit)),'rows');
        aMidbMid   = ismember(chComb,table2array(combinations(chSplit+1:2*chSplit,chSplit+1:2*chSplit)),'rows');
        aMidbDeep  = ismember(chComb,table2array(combinations(chSplit+1:2*chSplit,2*chSplit+1:maxChB)),'rows');

        aDeepbSuper = ismember(chComb,table2array(combinations(2*chSplit+1:maxChA,1:chSplit)),'rows');
        aDeepbMid   = ismember(chComb,table2array(combinations(2*chSplit+1:maxChA,chSplit+1:2*chSplit)),'rows');
        aDeepbDeep  = ismember(chComb,table2array(combinations(2*chSplit+1:maxChA,2*chSplit+1:maxChB)),'rows');        

        for iRef = 1:4
            switch iRef
                case 1
                    pacRef = pacVals(fileNum,iDate).aAmpBPhase;
                case 2
                    pacRef = pacVals(fileNum,iDate).aPhaseBAmp;
                case 3
                    pacRef = pacVals(fileNum,iDate).a2a; 
                    maxCh  = maxChA;
                    try
                        combSingleCh = allPACVars{fileNum,iDate}.combAllChA2A;
                    catch
                        combSingleCh = chComb; 
                    end

                case 4
                    pacRef = pacVals(fileNum,iDate).b2b;
                    maxCh  = maxChB;
                    try
                    combSingleCh = allPACVars{fileNum,iDate}.combAllChB2B;
                    catch
                        combSingleCh = chComb;
                    end

            end
            if iRef<3
                btwProbesCompMI(:,1) = reshape(squeeze(median(pacRef(:,:,:,aSuperbSuper),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,2) = reshape(squeeze(median(pacRef(:,:,:,aSuperbMid),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,3) = reshape(squeeze(median(pacRef(:,:,:,aSuperbDeep),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,4) = reshape(squeeze(median(pacRef(:,:,:,aMidbSuper),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,5) = reshape(squeeze(median(pacRef(:,:,:,aMidbMid),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,6) = reshape(squeeze(median(pacRef(:,:,:,aMidbDeep),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,7) = reshape(squeeze(median(pacRef(:,:,:,aDeepbSuper),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,8) = reshape(squeeze(median(pacRef(:,:,:,aDeepbMid),[1 4],'omitnan')),[nHigh*nLow 1]);
                btwProbesCompMI(:,9) = reshape(squeeze(median(pacRef(:,:,:,aDeepbDeep),[1 4],'omitnan')),[nHigh*nLow 1]);
            
            else
                superSuper = ismember(combSingleCh,table2array(combinations(1:chSplit,1:chSplit)),'rows');
                superMid   = ismember(combSingleCh,table2array(combinations(1:chSplit,chSplit+1:2*chSplit)),'rows');
                superDeep  = ismember(combSingleCh,table2array(combinations(1:chSplit,2*chSplit+1:maxCh)),'rows');

                midSuper = ismember(combSingleCh,table2array(combinations(chSplit+1:2*chSplit,1:chSplit)),'rows');
                midMid   = ismember(combSingleCh,table2array(combinations(chSplit+1:2*chSplit,chSplit+1:2*chSplit)),'rows');
                midDeep  = ismember(combSingleCh,table2array(combinations(chSplit+1:chSplit,2*chSplit+1:maxCh)),'rows');

                deepSuper = ismember(combSingleCh,table2array(combinations(2*chSplit+1:maxCh,1:chSplit)),'rows');
                deepMid   = ismember(combSingleCh,table2array(combinations(2*chSplit+1:maxCh,chSplit+1:2*chSplit)),'rows');
                deepDeep  = ismember(combSingleCh,table2array(combinations(2*chSplit+1:maxCh,2*chSplit+1:maxCh)),'rows');
                
                withinProbeComp(:,1) = reshape(squeeze(median(pacRef(:,:,:,superSuper),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,2) = reshape(squeeze(median(pacRef(:,:,:,superMid),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,3) = reshape(squeeze(median(pacRef(:,:,:,superDeep),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,4) = reshape(squeeze(median(pacRef(:,:,:,midSuper),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,5) = reshape(squeeze(median(pacRef(:,:,:,midMid),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,6) = reshape(squeeze(median(pacRef(:,:,:,midDeep),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,7) = reshape(squeeze(median(pacRef(:,:,:,deepSuper),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,8) = reshape(squeeze(median(pacRef(:,:,:,deepMid),[1 4],'omitnan')),[nHigh*nLow 1]);
                withinProbeComp(:,9) = reshape(squeeze(median(pacRef(:,:,:,deepDeep),[1 4],'omitnan')),[nHigh*nLow 1]);

            end
            switch iRef
                case 1
                    aAmpBPhaseLaminar(iDate,fileNum,:,:) = btwProbesCompMI; 
                case 2
                    aPhaseBAmpLaminar(iDate,fileNum,:,:) = btwProbesCompMI;
                case 3
                    a2aLaminar(iDate,fileNum,:,:) = withinProbeComp;
                case 4
                    b2bLaminar(iDate,fileNum,:,:) = withinProbeComp;
            end
        end
    end
end

%%
layerCombs = table2array(combinations(['S';'M';'D'],['S';'M';'D']));
aAmpBPhaseLaminarT = reshape(aAmpBPhaseLaminar,[matSize(1)*matSize(2) matSize(3) size(b2bLaminar,4)]); aAmpBPhaseLaminarT(zeroVals,:,:)=[];
aPhaseBAmpLaminarT = reshape(aPhaseBAmpLaminar,[matSize(1)*matSize(2) matSize(3) size(b2bLaminar,4)]);aPhaseBAmpLaminarT(zeroVals,:,:) =[];
a2aLaminarT        = reshape(a2aLaminar,[matSize(1)*matSize(2) matSize(3) size(b2bLaminar,4)]);a2aLaminarT(zeroVals,:,:) =[];
b2bLaminarT        = reshape(b2bLaminar,[matSize(1)*matSize(2) matSize(3) size(b2bLaminar,4)]);b2bLaminarT(zeroVals,:,:) = [];

figure; subplot(121); boxplot(squeeze(median(aAmpBPhaseLaminarT,2,'omitnan')),layerCombs); box off;
subplot(122);boxplot(squeeze(median(aPhaseBAmpLaminarT,2,'omitnan')),layerCombs); box off;

aAmpBPhaseThetaHighG = squeeze(median(aAmpBPhaseLaminarT(:,thetaHighGammaIdx,:),2,'omitnan')); 
aAmpBPhaseThetaLowG  = squeeze(median(aAmpBPhaseLaminarT(:,thetaLowGammaIdx,:),2,'omitnan'));

aPhaseBAmpThetaHighG = squeeze(median(aPhaseBAmpLaminarT(:,thetaHighGammaIdx,:),2,'omitnan')); 
aPhaseBAmpThetaLowG  = squeeze(median(aPhaseBAmpLaminarT(:,thetaLowGammaIdx,:),2,'omitnan'));

a2aLaminarThetaHighG = squeeze(median(a2aLaminarT(:,thetaHighGammaIdx,:),2,'omitnan')); 
a2aLaminarThetaLowG  = squeeze(median(a2aLaminarT(:,thetaLowGammaIdx,:),2,'omitnan'));

b2bLaminarThetaHighG = squeeze(median(b2bLaminarT(:,thetaHighGammaIdx,:),2,'omitnan')); 
b2bLaminarThetaLowG  = squeeze(median(b2bLaminarT(:,thetaLowGammaIdx,:),2,'omitnan'));

save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\allPACVars.mat'],...
    'aAmpBPhaseLaminarT','aPhaseBAmpLaminarT','a2aLaminarT','b2bLaminarT',...
    'aAmpBPhaseThetaHighG','aAmpBPhaseThetaLowG','aPhaseBAmpThetaHighG','aPhaseBAmpThetaLowG',...
    'a2aLaminarThetaHighG','a2aLaminarThetaLowG','b2bLaminarThetaHighG','b2bLaminarThetaLowG','-append');

%%
figure; subplot(121);boxplot([median(aAmpBPhaseThetaHighG(:,1:3),1,'omitnan')' median(aAmpBPhaseThetaHighG(:,4:6),1,'omitnan')' ...
    median(aAmpBPhaseThetaHighG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-high gamma')
 subplot(122);boxplot([median(aAmpBPhaseThetaLowG(:,1:3),1,'omitnan')' median(aAmpBPhaseThetaLowG(:,4:6),1,'omitnan')' ...
    median(aAmpBPhaseThetaLowG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-low gamma');

 %%
 figure; subplot(121);boxplot([median(aPhaseBAmpThetaHighG(:,1:3),1,'omitnan')' median(aPhaseBAmpThetaHighG(:,4:6),1,'omitnan')' ...
    median(aPhaseBAmpThetaHighG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-high gamma')
 subplot(122);boxplot([median(aPhaseBAmpThetaLowG(:,1:3),1,'omitnan')' median(aPhaseBAmpThetaLowG(:,4:6),1,'omitnan')' ...
    median(aPhaseBAmpThetaLowG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-low gamma');

 %%
 figure; subplot(121);boxplot([median(a2aLaminarThetaHighG(:,1:3),1,'omitnan')' median(a2aLaminarThetaHighG(:,4:6),1,'omitnan')' ...
    median(a2aLaminarThetaHighG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-high gamma')
 subplot(122);boxplot([median(a2aLaminarThetaLowG(:,1:3),1,'omitnan')' median(a2aLaminarThetaLowG(:,4:6),1,'omitnan')' ...
    median(a2aLaminarThetaLowG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-low gamma');

 %%
 figure; subplot(121);boxplot([median(b2bLaminarThetaHighG(:,1:3),1,'omitnan')' median(b2bLaminarThetaHighG(:,4:6),1,'omitnan')' ...
    median(b2bLaminarThetaHighG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-high gamma')
 subplot(122);boxplot([median(b2bLaminarThetaLowG(:,1:3),1,'omitnan')' median(b2bLaminarThetaLowG(:,4:6),1,'omitnan')' ...
    median(b2bLaminarThetaLowG(:,7:9),1,'omitnan')'],{'Superficial','Middle','Deep'});box off; title('Theta-low gamma');


 %%
 figure;
 subplot(121);boxplot([squeeze(median(aPhaseBAmpLaminar(:,1:3,:),[1 2],'omitnan')) squeeze(median(aPhaseBAmpLaminar(:,4:6,:),[1 2],'omitnan')) ...
    squeeze(median(aPhaseBAmpLaminar(:,7:9,:),[1 2],'omitnan'))],{'Superficial','Middle','Deep'});box off; title('Reference: Probe A Phase');

subplot(122);boxplot([squeeze(median(aPhaseBAmpLaminar(:,[1 4 7],:),[1 2],'omitnan')) squeeze(median(aPhaseBAmpLaminar(:,[2 5 8],:),[1 2],'omitnan')) ...
    squeeze(median(aPhaseBAmpLaminar(:,[3 6 9],:),[1 2],'omitnan'))],{'Superficial','Middle','Deep'});box off; title('Reference: Probe B Amplitude');
 
%%
 figure;
 subplot(121);boxplot([squeeze(median(aPhaseBAmpLaminar(:,1:3,:),[1 2],'omitnan')) squeeze(median(aPhaseBAmpLaminar(:,4:6,:),[1 2],'omitnan')) ...
    squeeze(median(aPhaseBAmpLaminar(:,7:9,:),[1 2],'omitnan'))],{'Superficial','Middle','Deep'});box off; title('Reference: Probe A Phase');

subplot(122);boxplot([squeeze(median(aPhaseBAmpLaminar(:,[1 4 7],:),[1 2],'omitnan')) squeeze(median(aPhaseBAmpLaminar(:,[2 5 8],:),[1 2],'omitnan')) ...
    squeeze(median(aPhaseBAmpLaminar(:,[3 6 9],:),[1 2],'omitnan'))],{'Superficial','Middle','Deep'});box off; title('Reference: Probe B Amplitude');


%%
 figure;
 subplot(121);boxplot([squeeze(median(aAmpBPhaseLaminar(:,1:3,:),[1 2],'omitnan')) squeeze(median(aAmpBPhaseLaminar(:,4:6,:),[1 2],'omitnan')) ...
    squeeze(median(aAmpBPhaseLaminar(:,7:9,:),[1 2],'omitnan'))],{'Superficial','Middle','Deep'});box off; title('Reference: Probe A Phase');

subplot(122);boxplot([squeeze(median(aAmpBPhaseLaminar(:,[1 4 7],:),[1 2],'omitnan')) squeeze(median(aAmpBPhaseLaminar(:,[2 5 8],:),[1 2],'omitnan')) ...
    squeeze(median(aAmpBPhaseLaminar(:,[3 6 9],:),[1 2],'omitnan'))],{'Superficial','Middle','Deep'});box off; title('Reference: Probe B Amplitude');


%%

aAmpBPhaseThetaLowG = squeeze(median(aAmpBPhaseLaminar(thetaLowGammaIdx,:,:),1,'omitnan'));
aAmpBPhaseThetaHighG = squeeze(median(aAmpBPhaseLaminar(thetaHighGammaIdx,:,:),1,'omitnan'));

aPhaseBAmpThetaLowG = squeeze(median(aPhaseBAmpLaminar(thetaLowGammaIdx,:,:),1,'omitnan'));
aPhaseBAmpThetaHighG = squeeze(median(aPhaseBAmpLaminar(thetaHighGammaIdx,:,:),1,'omitnan'));


% Change reference between amplitude and phase
ampRefThetaLowG = [median([aAmpBPhaseThetaLowG(1:3,:);aPhaseBAmpThetaLowG([1 4 7],:)],1,'omitnan')'  ...
    median([aAmpBPhaseThetaLowG(4:6,:);aPhaseBAmpThetaLowG([2 5 8],:)],1,'omitnan')' ...
    median([aAmpBPhaseThetaLowG(7:9,:);aPhaseBAmpThetaLowG([3 6 9],:)],1,'omitnan')'];

ampRefThetaHighG = [median([aAmpBPhaseThetaHighG(1:3,:);aPhaseBAmpThetaHighG([1 4 7],:)],1,'omitnan')'  ...
    median([aAmpBPhaseThetaHighG(4:6,:);aPhaseBAmpThetaHighG([2 5 8],:)],1,'omitnan')' ...
    median([aAmpBPhaseThetaHighG(7:9,:);aPhaseBAmpThetaHighG([3 6 9],:)],1,'omitnan')'];

phaseRefThetaLowG = [median([aAmpBPhaseThetaLowG([1 4 7],:);aPhaseBAmpThetaLowG(1:3,:)],1,'omitnan')'  ...
    median([aAmpBPhaseThetaLowG([2 5 8],:);aPhaseBAmpThetaLowG(4:6,:)],1,'omitnan')' ...
    median([aAmpBPhaseThetaLowG([3 6 9],:);aPhaseBAmpThetaLowG(7:9,:)],1,'omitnan')'];

phaseRefThetaHighG = [median([aAmpBPhaseThetaHighG([1 4 7],:);aPhaseBAmpThetaHighG(1:3,:)],1,'omitnan')'  ...
    median([aAmpBPhaseThetaHighG([2 5 8],:);aPhaseBAmpThetaHighG(4:6,:)],1,'omitnan')' ...
    median([aAmpBPhaseThetaHighG([3 6 9],:);aPhaseBAmpThetaHighG(7:9,:)],1,'omitnan')'];
%
for iRef = 1:2
    switch iRef
        case 1
            thetaLG = aAmpBPhaseThetaLowG;
            thetaHG = aAmpBPhaseThetaHighG;
            plTitle = 'A-amplitude B-Phase';
        case 2
            thetaLG = aPhaseBAmpThetaLowG;
            thetaHG = aPhaseBAmpThetaHighG;
            plTitle = 'A-phase B-Amplitude';
    end
    figure;
    subplot(121); boxplot(thetaLG',layerCombs); title(['Theta-low Gamma: ' plTitle]); box off; ylim([2e-4 6e-4]);
    subplot(122); boxplot(thetaHG',layerCombs); title(['Theta-High Gamma: ' plTitle]); box off;ylim([3.3e-4 6.5e-4]);

    figure; 
    subplot(121); boxplot([median(thetaLG(1:3,:),1,'omitnan')' median(thetaLG(4:6,:),1,'omitnan')' ...
        median(thetaLG(7:9,:),1,'omitnan')'],{'Superficial','Middle','Deep'});
    title(['Theta-Low Gamma: ' plTitle]);box off; ylim([2e-4 6e-4]);
     subplot(122); boxplot([median(thetaHG(1:3,:),1,'omitnan')' median(thetaHG(4:6,:),1,'omitnan')' ...
        median(thetaHG(7:9,:),1,'omitnan')'],{'Superficial','Middle','Deep'});
     title(['Theta-High Gamma: ' plTitle]); box off;ylim([3.3e-4 6e-4]);
   
end

%% Figures
figure; subplot(221); boxplot(ampRefThetaLowG,{'Superficial','Middle','Deep'}); box off; axis square; title('Theta-Low Gamma - Amplitude Ref'); ylim([2e-4 5.5e-4]);
subplot(222); boxplot(ampRefThetaHighG,{'Superficial','Middle','Deep'}); box off; axis square; title('Theta-High Gamma - Amplitude Ref'); ylim([3.5e-4 5.8e-4]);
subplot(223); boxplot(phaseRefThetaLowG,{'Superficial','Middle','Deep'}); box off; axis square; title('Theta-Low Gamma - Phase Ref'); ylim([2e-4 5.5e-4]);
subplot(224); boxplot(phaseRefThetaHighG,{'Superficial','Middle','Deep'}); box off; axis square; title('Theta-High Gamma - Phase Ref'); ylim([3.5e-4 5.8e-4]);

%% Laminar results - need to move to combinedAnalysisEphysDualProbe
% pairClass = allVars.pairClass; 
singleChRow = cellfun(@(x) isscalar(x),allVars.intraCorrBR(:,1));

pairClass(singleChRow,:) = [];

% Time series
superAll      = allVars.medPairCorrSuperR; 
midAll        = allVars.medPairCorrMidR;
deepAll       = allVars.medPairCorrDeepR; 
superMidPair  = (allVars.medPairASuperBMidR + allVars.medPairAMidBSuperR)./2;
superDeepPair = (allVars.medPairASuperBDeepR + allVars.medPairADeepBSuperR)./2 ;
midDeepPair   = (allVars.medPairAMidBDeepR + allVars.medPairADeepBMidR)./2 ; 

% Power
superPow         = allVars.envelopePairCorrSuperR; 
midPow           = allVars.envelopePairCorrMidR; 
deepPow          = allVars.envelopePairCorrDeepR; 
superMidPairPow  = (allVars.envelopeASuperBMidR + allVars.envelopeAMidBSuperR)./2 ;
superDeepPairPow = (allVars.envelopeADeepBSuperR + allVars.envelopeASuperBDeepR)./2 ; 
midDeepPairPow   = (allVars.envelopeAMidBDeepR + allVars.envelopeADeepBMidR)./2 ;

% Infraslow
superInfra         = allVars.infraPairCorrSuperR; 
midInfra           = allVars.infraPairCorrMidR; 
deepInfra          = allVars.infraPairCorrDeepR; 
superMidPairInfra  = (allVars.infraASuperBMidR + allVars. infraAMidBSuperR)./2 ;
superDeepPairInfra = (allVars.infraASuperBDeepR + allVars.infraADeepBSuperR)./2 ; 
midDeepPairInfra   = (allVars.infraAMidBDeepR + allVars.infraADeepBMidR)./2 ;

connValsR = allVars.connValsR; connValsR(singleChRow) = [];
distValsR = allVars.distValsR; distValsR(singleChRow) = [];

smLoc = sum(pairClass=='SM',2)==2;
ssLoc = sum(pairClass=='SS',2)==2;
mmLoc = sum(pairClass=='MM',2)==2;

%% Identify sites where termination is in 3b or 4
termination3b = false(max(cell2mat(cellfun(@length,datFileNumAll,'un',0))),size(allDates,1)); 
terminationArea4 = false(max(cell2mat(cellfun(@length,datFileNumAll,'un',0))),size(allDates,1));
for iDate = 1:size(allDates,1)
    datFileNum = datFileNumAll{iDate,1};
    for iRun = 1:length(datFileNum)
        fileNum = datFileNum(iRun);

        locs{fileNum,iDate}  = {corticalAreaProbeA{iDate}{fileNum} corticalAreaProbeB{iDate}{fileNum} };

        % Termination in 3b (1-->3b, 2-->3b)
        if strcmp(corticalAreaProbeA{iDate}(fileNum),'3b') || strcmp(corticalAreaProbeB{iDate}(fileNum),'3b')
            if strcmp(corticalAreaProbeA{iDate}(fileNum),'3b')
                if strcmp(corticalAreaProbeB{iDate}(fileNum),'1') || strcmp(corticalAreaProbeB{iDate}(fileNum),'2') || strcmp(corticalAreaProbeB{iDate}(fileNum),'3a')
                    termination3b(fileNum,iDate) = true;
                end

            else
                if strcmp(corticalAreaProbeA{iDate}(fileNum),'1') || strcmp(corticalAreaProbeA{iDate}(fileNum),'2')|| strcmp(corticalAreaProbeA{iDate}(fileNum),'3a')
                    termination3b(fileNum,iDate) = true;           
                end
            end
        end

        % Pairs of 6 and 4
        if strcmp(corticalAreaProbeA{iDate}(fileNum),'4') || strcmp(corticalAreaProbeB{iDate}(fileNum),'4')
            if strcmp(corticalAreaProbeA{iDate}(fileNum),'4')
                if strcmp(corticalAreaProbeB{iDate}(fileNum),'6')
                    terminationArea4(fileNum,iDate) = true;
                end
            else
                if strcmp(corticalAreaProbeA{iDate}(fileNum),'6')
                    terminationArea4(fileNum,iDate) = true;
                end
            end
        end
    end
end
if strcmp(monkeyName,'CharlieSheen')
    loc3bInfo = zeros(size(termination3b));
    loc3bInfo(8:11,1) = 1; loc3bInfo(5:6,3)= 1; 
    loc3bInfo(1:5,4) = 2; loc3bInfo([3 5],5)= 2; 
    loc3bInfo([1:3 5],6)= 1; loc3bInfo(15:16,6)= 2; 
else
    loc3bInfo = termination3b;
end

loc3bInfo = reshape(loc3bInfo',[numel(loc3bInfo) 1]);
loc3bInfo(allVars.removeDataIdx) = []; loc3bInfo(singleChRow) = [];

termination3b = reshape(termination3b',[numel(termination3b) 1]); 
termination3b(allVars.removeDataIdx) = []; termination3b(singleChRow) = [];

termination3bARef = zeros(size(termination3b)); 
termination3bARef(loc3bInfo==1) = 1;


termination3bBRef = zeros(size(termination3b));
termination3bBRef(loc3bInfo==2) = 1; 

terminationArea4 = reshape(terminationArea4',[numel(terminationArea4) 1]); 
terminationArea4(allVars.removeDataIdx) = []; terminationArea4(singleChRow) = [];

save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],...
    'termination3bARef','termination3bBRef','-append');

%% Show the distribution of correlations with respect to a reference compartment
flagVal = mmLoc;%true(size(connValsR)); 
figure; idx = 1;
for iRef = 1:3
    switch iRef
        case 1
            superVal = superAll;
            midVal   = midAll;
            deepVal  = deepAll;
            superMid = superMidPair;
            superDeep = superDeepPair; 
            midDeep   = midDeepPair; 
            typeLabel = 'Time series';

        case 2
            superVal = superPow ;
            midVal   = midPow ;
            deepVal  = deepPow ;
            superMid = superMidPairPow;
            superDeep = superDeepPairPow ; 
            midDeep   = midDeepPairPow ; 
            typeLabel = 'Power';

        case 3
            superVal = superInfra;
            midVal   = midInfra;
            deepVal  = deepInfra;
            superMid = superMidPairInfra;
            superDeep = superDeepPairInfra; 
            midDeep   = midDeepPairInfra; 
            typeLabel = 'Infraslow';
    end

    [superAllCorr(iRef,:),pValSuperAll(iRef,:)] = partialcorr(superVal(flagVal,:),connValsR(flagVal),distValsR(flagVal));
    [midAllCorr(iRef,:),pValMidAll(iRef,:)]     = partialcorr(midVal(flagVal,:),connValsR(flagVal),distValsR(flagVal));
    [deepAllCorr(iRef,:),pValDeepAll(iRef,:)]   = partialcorr(deepVal(flagVal,:),connValsR(flagVal),distValsR(flagVal));

    [superMidCorr(iRef,:),pValsSuperMid(iRef,:)]   = partialcorr(superMid(flagVal,:),connValsR(flagVal),distValsR(flagVal));
    [superDeepCorr(iRef,:),pValsSuperDeep(iRef,:)] = partialcorr(superDeep(flagVal,:),connValsR(flagVal),distValsR(flagVal));
    [midDeepCorr(iRef,:),pValsMidDeep(iRef,:)]     = partialcorr(midDeep(flagVal,:),connValsR(flagVal),distValsR(flagVal));

    [superRefCorr(iRef,:),pSuperRef(iRef,:)] = partialcorr([superVal(flagVal,:); superMid(flagVal,:); superDeep(flagVal,:)],...
        repmat(connValsR(flagVal),[3 1]),repmat(distValsR(flagVal),[3 1]));

    [midRefCorr(iRef,:),pMidRefCorr(iRef,:)] = partialcorr([midVal(flagVal,:); superMid(flagVal,:); midDeep(flagVal,:)],...
        repmat(connValsR(flagVal),[3 1]),repmat(distValsR(flagVal),[3 1]));

     [deepRefCorr(iRef,:),pDeepRefCorr(iRef,:)] = partialcorr([deepVal(flagVal,:); midDeep(flagVal,:); superDeep(flagVal,:)],...
        repmat(connValsR(flagVal),[3 1]),repmat(distValsR(flagVal),[3 1]));


     subplot(3,3,idx); boxplot([superVal(flagVal,:); superMid(flagVal,:); superDeep(flagVal,:)]);
     xticklabels(bandLabels); title(['Superficial - ' typeLabel]);box off;yticks(-0.4:0.2:1); axis square;
     if iRef~=2; ylim([-0.4 1]);  yticks(-0.4:0.2:1); else; yticks([-0.2 0.6]);  yticks(-0.2:0.1:0.6);end

     subplot(3,3,idx+1); boxplot([midVal(flagVal,:); superMid(flagVal,:); midDeep(flagVal,:)]);
     xticklabels(bandLabels); title(['Middle - ' typeLabel]);box off;yticks(-0.4:0.2:1);axis square;
     if iRef~=2; ylim([-0.4 1]);  yticks(-0.4:0.2:1); else; yticks([-0.2 0.6]); yticks(-0.2:0.1:0.6); end

     subplot(3,3,idx+2); boxplot([deepVal(flagVal,:); midDeep(flagVal,:); superDeep(flagVal,:)]);
     xticklabels(bandLabels); title(['Deep - ' typeLabel]);box off;yticks(-0.4:0.2:1);axis square;
     if iRef~=2; ylim([-0.4 1]);  yticks(-0.4:0.2:1); else; yticks([-0.2 0.6]); yticks(-0.2:0.1:0.6); end

     idx = idx+3;

end

%%
figure; 
subplot(231); imagesc(1:5,1:3,superAllCorr); title('S/S');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.2 0.5]); colorbar;
subplot(232); imagesc(1:5,1:3,superMidCorr); title('S/M');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.2 0.5]);colorbar;
subplot(233); imagesc(1:5,1:3,superDeepCorr); title('S/D');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.2 0.5]);colorbar;
subplot(234); imagesc(1:5,1:3,midAllCorr); title('M/M');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.2 0.5]);colorbar;
subplot(235); imagesc(1:5,1:3,midDeepCorr); title('M/D');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.2 0.5]);colorbar;
subplot(236); imagesc(1:5,1:3,deepAllCorr); title('D/D');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.2 0.5]);colorbar;

%%
figure;
subplot(131); imagesc(1:5,1:3,superRefCorr); title('Superficial');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.1 0.5]); colorbar;
subplot(132); imagesc(1:5,1:3,midRefCorr); title('Middle');axis square; xticklabels( bandLabels); yticklabels(timeLabels); clim([-0.1 0.5]); colorbar;
subplot(133); imagesc(1:5,1:3,deepRefCorr); title('Deep');axis square; xticklabels(bandLabels); yticklabels(timeLabels); clim([-0.1 0.5]); colorbar;



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

function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off;
coeff = polyfit(xVal,yVal,1);
xFit  = linspace(min(xVal),max(xVal),1000);
yFit  = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
% coeff = fit(xVal,double(yVal),'poly1','Robust','LAR');
% xFit  = linspace(min(xVal),max(xVal),1000);
% yFit  = coeff.p1*xFit + coeff.p2; mdl = fitlm(xVal,yVal,'RobustOpts','on');
plot(xFit,yFit,'-k','LineWidth',1);
if nargin<3
    textLocX  = max(xVal)-0.2*max(xVal);
    textLocY1 = max(yVal)-0.2*max(yVal);
    textLocY2 = max(yVal)-0.3*max(yVal);
end
text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end

function [fVals,psd1]  = getFFTAllCh(signal,fs)

l = size(signal,1);
fVals = fs*(0:round(l/2))/l;
psd2  =  abs(fft(signal,[],2)/l);
psd1  = psd2(1:round(l/2)+1,:); 
psd1(2:end-1,:) = 2*psd1(2:end-1,:); 

end

function [PSD1, freq] = calculatePSD(sig,params)
if size(sig,2)==1
    [PSD,freq] = mtspectrumc(buffer(sig,2048*2),params);
    PSD1 = median(PSD,2);
else
    [PSD1,freq] = mtspectrumc(sig,params);
end
end