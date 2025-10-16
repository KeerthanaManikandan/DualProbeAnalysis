function [allProbeData,allBadTimes,badElecA,badElecB,estChInCortexA,estChInCortexB] = getAllDualProbeData(monkeyName,hemisphere,allDates, ...
    datFileNameAll, datFileNumAll,serverPath,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag )
% This function stores/retrieves LFP data, identifies bad channels,
% identifies bad time segments, and transition channels for all recording
% pairs. 
% Keerthana Manikandan- October 15, 2025

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

% Determine bad time segments and channels from the two probes
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\badChTimeVars.mat'],'file') 
    tic;
    [allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,allProbeData,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag);
    toc;
    disp(['Identified bad channels and time segments for ' monkeyName]);
    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\badChTimeVars.mat'],'allBadTimes','badElecA','badElecB');
else
    load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\badChTimeVars.mat']);
    disp(['Retrieving bad channels and time segments for ' monkeyName]);
end

% Get the transition channels for the recordings...
% Eliminate bad channels and obtain the transition channels by
% computing slope of the marginals (obtained by averaging the intra-probe
% correlograms) - Mean is used here for sensitivity to outliers which is
% needed for picking the transition channel
clear estChInCortexA estChInCortexB

 disp(['Identifying transition channels for ' monkeyName]);
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

end