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

fs             = 1e3; % Sampling frequency 
gammaBand      = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); 
alphaBand      = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); 
betaBand       = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
hemisphere     = 'Left';
saveFigureFlag = 0;


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

[allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate,patchInfo,pairClass] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']);

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

            [probe1Ch,probe2Ch,eegCh,scalpEEGCh] = ...
                saveLFPDualProbe(monkeyName,expDate,fileNum,datFileName,saveFolder,serverPath,fs);

            % Save the LFP
            if ~exist('saveFolder','dir'); [~,~] = mkdir(saveFolder); end
            disp('Storing data... ' );

            if exist('eegCh','var') && exist('scalpEEGCh','var')
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch','eegCh','scalpEEGCh');

            elseif exist('eegCh','var') && ~exist('scalpEEGCh','var')
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch','eegCh');

            elseif  ~(exist('eegCh','var') || exist('scalpEEGCh','var'))
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probe1Ch','probe2Ch');
            end
            
            allProbeData{fileNum,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);

        else
            % Retrieve LFP Data
            disp(['Retrieving data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            allProbeData{fileNum,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
        
            % probe1{fileNum,iDate} = probe1Ch;
            % probe2{fileNum,iDate} = probe2Ch;
            % 
            % if exist('eegCh','var') % Get EEG
            %     eeg{fileNum,iDate} = single(eegCh);
            % else
            %     eeg{fileNum,iDate} = [];
            % end
            % 
            % if exist('scalpEEGCh','var') % Get Scalp EEG
            %     scalpEEG{fileNum,iDate} = single(scalpEEGCh);
            % else
            %     scalpEEG{fileNum,iDate} = [];
            % end
        end
    end
end

clear probe1Ch probe2Ch eegCh
clc; disp('Data Stored/Retrieved');

%% Determine and remove bad time segments and channels from the two probes
% Determine bad time segments and bad channels
tic;
[allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,allProbeData,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);
toc;

% Remove bad channels and times
% for iDate =1:size(allDates,1) % all experiment dates
%     clear expDate datFileNum
%     expDate    = allDates(iDate,:);
%     datFileNum = datFileNumAll{iDate,1};
% 
%     for iRun = 1:length(datFileNum) % all recordings for a particular experiment date
%         if strcmp(expDate,'09_19_2022') && iRun == 4; continue; end
%         if strcmp(expDate,'02_07_2023') && iRun == 2; continue; end
%         if strcmp(expDate,'04_11_2023') && iRun ==10; continue; end
%         fileNum = datFileNum(iRun);
% 
%         if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP
% 
%         % Remove extra channels if present
%         if size(probe1{fileNum,iDate},2) == 33; probe1{fileNum,iDate}(:,1) = []; end
%         if size(probe2{fileNum,iDate},2) == 33; probe2{fileNum,iDate}(:,1) = []; end
% 
%         % Remove bad channels from both probes
%         if ~isempty(badElecA{fileNum,iDate})
%             probe1{fileNum,iDate}(:,badElecA{fileNum,iDate}) = [];
%         end
% 
%         if ~isempty(badElecB{fileNum,iDate})
%             probe2{fileNum,iDate}(:,badElecB{fileNum,iDate}) = [];
%         end
% 
%         % Remove the bad time segments
%         if ~isempty(allBadTimes{fileNum,iDate})
%             probe1{fileNum,iDate}(allBadTimes{fileNum,iDate},:) = [];
%             probe2{fileNum,iDate}(allBadTimes{fileNum,iDate},:) = [];
%         end
%     end
% end
clc; disp(['Identified and removed bad channels, time segments for ' monkeyName]);

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
        if strcmp(expDate,'09_19_2022') && iRun == 4; continue; end
        if strcmp(expDate,'02_07_2023') && iRun == 2; continue; end
        if strcmp(expDate,'04_11_2023') && iRun ==10; continue; end
        fileNum = datFileNum(iRun);

        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate}); continue; end % Check if probes do not have any LFP

        % Get the envelope (Wideband power) for each electrode 
        clear probeA probeB marginalVal fx
        probeA = envelope(probe1{fileNum,iDate},5);
        probeB = envelope(probe2{fileNum,iDate},5);

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

                % if ~exist(fullfile(dataDir,'pairCorr.png'),'file') % Within probe correlations
                %     tiledlayout(1,3,'TileSpacing','Compact');
                %     nexttile([1 2]); imagesc(imgaussfilt(corr(probeTemp,'Rows','complete'),1));
                %     colormap jet; colorbar; axis image tight;
                %     xticks(1:2:size(probeTemp,2));  yticks(1:2:size(probeTemp,2));
                % 
                %     nexttile; plot(marginalVal,1:length(marginalVal));set(gca,'YDir','reverse');
                %     axis padded;xlim([0 round((max(marginalVal)+0.1).*10)/10]);
                %     xticks(0:0.2:round((max(marginalVal)+0.1).*10)/10);
                % 
                %     sgtitle(strrep(['Within probe correlations for: ' monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
                %     f = gcf; exportgraphics(f,[dataDir '\pairCorr.png'],'Resolution',300); close gcf;
                % end

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
clc; disp(['Determined transition channels for ' monkeyName]);

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