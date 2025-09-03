function [allBadTimes,badElecA,badElecB] = getBadTimesAndChannels(monkeyName, allDates, datFileNumAll,allProbeData,chInCortexProbeA,chInCortexProbeB,probeLabelA,probeLabelB,saveFigureFlag)
% This function finds and removes the bad time segments and bad channels
% and also outputs the spectrograms (mean and median across channels)
% Update log:
% Date: March 1, 2024: KM to modify the code for clarity and consistency

if ~exist('saveFigureFlag','var') || isempty(saveFigureFlag); saveFigureFlag = 0; end
fs              = 1e3; % Sampling frequency
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

hemisphere      = 'Left';
chOutCortex     = 1:3;
chDeep          = 30:32;

gammaBand      = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand      = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters

for iDate = 1:size(allDates,1) % All experiments
    clear expdate datFileNum datFileName
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];

    % Create folders
    if ~exist([saveFolder '\AvgSpectrograms'],'dir'); [~,~] = mkdir([saveFolder '\AvgSpectrograms']); end
    if ~exist([saveFolder '\IntraProbeCorr'],'dir');  [~,~] = mkdir([saveFolder '\IntraProbeCorr']);  end
    if ~exist([saveFolder '\InterProbeCorr'],'dir');  [~,~] = mkdir([saveFolder '\InterProbeCorr']);  end
    if ~exist([saveFolder '\Transition'],'dir');      [~,~] = mkdir([saveFolder '\Transition']);      end

    for iFile = 1:length(datFileNum)
        clear fileNum probe1 probe2 badTimesA badTimesB
        % if strcmp(expDate,'09_19_2022') && iFile == 4;  continue; end
        % if strcmp(expDate,'02_07_2023') && iFile == 2;  continue;  end
        % if strcmp(expDate,'04_11_2023') && iFile == 10; continue; end

        fileNum = datFileNum(iFile);
        probe1 = single(allProbeData{fileNum,iDate}.probe1Ch);
        probe2 = single(allProbeData{fileNum,iDate}.probe2Ch);

        if iDate>=6
            [bL,aL] = butter(3,([6 250]./(fs/2)),'bandpass');
            probe1 = single(filtfilt(bL,aL,double(probe1)));
            probe2 = single(filtfilt(bL,aL,double(probe2)));
        end

        disp(['Preprocessing: ' num2str(expDate) ' File: ' num2str(fileNum)]);

        [~,timeValsSpec1,~] = mtspecgramc(probe1,[5 2],params); % Checking spectrogram for bad recording/run 
        [~,timeValsSpec2,~] = mtspecgramc(probe2,[5 2],params);
        if isempty(timeValsSpec1)|| isempty(timeValsSpec2); continue; end   

        % Finding bad time segments
        for iProbe = 1:2
            clear probe channls spec timeVals powTimeBin badElecThresh badChVal badTimeThresh badTimes probeLabel

            if isempty(probe1) || isempty(probe2)
                badElecA{fileNum,iDate} = [];
                badElecB{fileNum,iDate} = [];

                continue;
            end
                        
            switch iProbe
                case 1 % Probe A
                    probe     = probe1;  
                    probeName = 'A';
                    if size(probe,2) == 33; probe(:,1) = []; end
                    channels  = 1:size(probe,2);

                    if ~isempty(probeLabelA{iDate,1})
                        probeLabel = probeLabelA{iDate,1}(fileNum,:);
                    else
                        probeLabel = [];
                    end                    

                case 2 % Probe B
                    probe     = probe2;
                    probeName = 'B';
                    channels  = 1:size(probe,2);

                    if ~isempty(probeLabelB{iDate,1})
                        probeLabel = probeLabelB{iDate,1}(fileNum,:);
                    else
                        probeLabel = [];
                    end

                    % if iDate>5
                    %     if size(probe,2) == 33; probe(:,1) = []; end
                    %     channels = [1:13 15:21 23:size(probe,2)];

                    % 
                    % elseif strcmp(expDate,'08_08_2022')
                    %     badTimesB = [];
                    %     probeLabel = [];
                    %     badElecB{fileNum,iDate} = [];
                    %     continue;
                    % else
                    %     channels = 1:size(probe,2);
                    %     if ~isempty(probeLabelB{iDate,1})
                    %         probeLabel = probeLabelB{iDate,1}(fileNum,:);
                    %     else
                    %         probeLabel = [];
                    %     end
                    % 
                    % end
            end                     

            % Remove bad channels
            badChVal = [];
            
            if size(probe,2)~=1
                
                % [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probe,[5 2],params);
                % freqIdx = freqValsSpec>=65 & freqValsSpec<=85;
                % powTimeBin = squeeze(sum(10.*log10(abs(spec(:,freqIdx,:))),2));
                % % [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probe,[5 2],params);
                % % freqIdx = freqValsSpec>=65 & freqValsSpec<=85;
                % % powTimeBin = squeeze(sum(spec(:,freqIdx,:),2));
                % 
                % if isempty(timeValsSpec); continue; end
                % 
                % if length(channels)~= 1       % Determine the minimum and maximum threshold to determine bad channels
                %     badElecThreshHigh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
                %     badElecThreshLow  = (median(powTimeBin,2)-5*mad(powTimeBin,1,2));
                % 
                %     badChVal = [badChVal channels(sum((powTimeBin>badElecThreshHigh),1) >= floor(0.75*size(powTimeBin,1)) |(sum((powTimeBin<badElecThreshLow),1) >= floor(0.75*size(powTimeBin,1))))];
                % end

                if ~isempty(probeLabel) % Check if previously known bad channels have been identified
                    if strcmp(probeLabel,'BD29')
                        if ~any(badChVal==8); badChVal = unique([badChVal 8]); end

                    elseif strcmp(probeLabel,'CDE1')
                        if ~any(badChVal==9); badChVal = unique([badChVal 9]); end

                    elseif strcmp(probeLabel,'D553')
                        if ~any(badChVal==5); badChVal = unique([badChVal 5]); end

                    elseif strcmp(probeLabel,'B25A')
                        if ~any(badChVal==27); badChVal = unique([badChVal 27]); end

                    elseif strcmp(probeLabel,'0763')
                        if ~(any(badChVal==2) && any(badChVal==4)); badChVal = unique([badChVal 2 4]); end
                    end
                end

                % Check if channels 14,22 are included in the bad channels
                % when nano2+stim front ends are used for recording.
                if iProbe==2 && ~(strcmp(expDate,'10_21_2024')|| strcmp(expDate,'04_30_2025'))
                    if sum(ismember(badChVal,[14 22]))~=2
                        badChVal = unique([badChVal 14 22]);
                    end
                end

                % Remove other channels that were identified from
                % spectrogram/raw power method
                if iProbe==2
                    if strcmp(expDate,'02_21_2023') && fileNum == 7
                        badChVal = unique([badChVal 29]);
                    elseif strcmp(expDate,'10_17_2022') && fileNum == 9
                        badChVal = unique([badChVal 27]);
                    end
                end

                % Remove bad channels from the probe
                if ~isempty(badChVal); channels(ismember(channels,badChVal)) = []; end

            else
                badChVal = [];
            end

            % Determine the bad time segments
            % clear probeBL
            % probeBL = single(filtfilt(bG,aG,double(probe(:,channels)))); % Filtering the gamma band
            %
            % % Bad time segment threshold bounds
            % if length(channels)~= 1 % Linear array
            %     badTimeThreshHigh = mean(probeBL(:,15:end),'all','omitnan') + 5*std(probeBL(:,15:end),[],[1,2]);
            %     badTimeThreshLow  = mean(probeBL(:,15:end),'all','omitnan') - 5*std(probeBL(:,15:end),[],[1,2]);
            % else % Single probe
            %     badTimeThreshHigh = mean(probeBL,'omitnan') + 5*std(probeBL,[]);
            %     badTimeThreshLow= mean(probeBL,'omitnan')- 5*std(probeBL,[]);
            % end
            %
            % % Determine threshold crossings...
            % badTimeIndOld = find(mean(probeBL(:,15:end),2,'omitnan')>badTimeThreshHigh | mean(probeBL(:,15:end),2,'omitnan')<badTimeThreshLow);
            %
            % badTimeInd = []; badTimes = [];
            %
            % if ~isempty(badTimeIndOld)
            %     % Taking 250 ms before and after each threshold crossing
            %     badTimeInd =[(badTimeIndOld-250)  (badTimeIndOld+250)];
            %
            %     for iL = 1:size(badTimeInd,1)
            %         badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
            %     end
            %
            %     badTimes = unique(badTimes);
            %     badTimes(badTimes>size(probeBL,1)) = [];
            %     badTimes(badTimes<=0) = [];
            % end
            clear probeBL badTimeThreshHigh badTimeThreshLow badTimeIndOld

            probeBL = envelope(probe(:,channels));
            chLim = 15; thresh = 4; timeBin = 50;

            badTimeThreshHigh = median(probeBL(:,chLim:end),'all','omitnan') +thresh*mad(probeBL(:,chLim:end),[],[1,2]);
            badTimeThreshLow  = median(probeBL(:,chLim:end),'all','omitnan') -thresh*mad(probeBL(:,chLim:end),[],[1,2]);
            badTimeIndOld = find(mean(probeBL(:,chLim:end),2,'omitnan')>=badTimeThreshHigh | mean(probeBL(:,chLim:end),2,'omitnan')<=badTimeThreshLow);
            badTimeInd = []; badTimes = []; clear badTimesCell

            badTimeInd   = [(badTimeIndOld-timeBin)  (badTimeIndOld+timeBin)];
            badTimesCell = arrayfun(@(x) badTimeInd(x,1):badTimeInd(x,2), 1:size(badTimeInd,1), 'UniformOutput', false);
            badTimes     = unique([badTimesCell{:}]);

            badTimes(badTimes>size(probeBL,1)) = [];
            badTimes(badTimes<=0)              = [];
            
            % Remove bad times as determined from the spectrogram, and that
            % did not get removed properly after bad time segment removal
            if iProbe==2 && strcmp(expDate,'11_01_2021') && fileNum ==10
                badTimes = unique([badTimes 210e3:310e3]);
                
            elseif iProbe == 2 && strcmp(expDate,'10_10_2022') && fileNum ==7
                badTimes = unique([badTimes 360e3:400e3]);
            end


            % [spec,timeValsSpec,~] = mtspecgramc(probe,[5 2],params);
            % meanS = mean(10.*log10(abs(spec(:,:,15:end))),3,'omitnan');
            % powMeanS = squeeze(sum(meanS,2));
            % 
            % badTimeThreshHigh   = (median(powMeanS,1)+3*mad(powMeanS,1,1));
            % badTimeThreshLow   = (median(powMeanS,1)-3*mad(powMeanS,1,1));
            % 
            % badTimeInd = floor(timeValsSpec(powMeanS>badTimeThreshHigh | powMeanS<badTimeThreshLow)*1e3);
            % 
            % badTimes = [];
            % if ~isempty(badTimeInd)
            %     jumpIdx = find(isoutlier([0 diff(badTimeInd)]));
            % 
            %     for iL = 1:length(jumpIdx)
            %         startIdx = badTimeInd(jumpIdx(iL));
            % 
            %         if iL == length(jumpIdx)
            %             endIdx = badTimeInd(end);
            %         else
            %             endIdx = badTimeInd(jumpIdx(iL+1)-1);
            %         end
            % 
            %         if startIdx == endIdx
            %             badTimes =[badTimes startIdx-50:startIdx+50];
            %         else
            %             badTimes = [badTimes startIdx:endIdx];
            %         end
            %     end
            %     badTimes = unique(badTimes);
            % end

            switch iProbe
                case 1
                    badTimesA               = badTimes;
                    badElecA{fileNum,iDate} = badChVal;
                case 2
                    badTimesB               = badTimes;
                    badElecB{fileNum,iDate} = badChVal;
            end


            % Average spectrogram before and after removal of bad time
            % segments...
            if ~exist([ saveFolder '\AvgSpectrograms\AvgSpecgram_Before_After_' num2str(iFile) '_Probe' probeName '.png'],'file') ||  saveFigureFlag
                figure('units','normalized','outerposition',[0 0 1 1]);

                subplot(2,1,1);
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probe(:,channels(channels>=15)),[5 2],params);
                meanS = mean(10.*log10(abs(spec)),3,'omitnan'); % Mean across channels 15-32
                imagesc(timeValsSpec,freqValsSpec,meanS'); colormap jet; shading interp; set(gca,'YDir','normal');
                xlabel('Time (s)'); ylabel('Frequency (Hz)'); clim([-30 30]); title('Spectrogram before removal of bad time segments'); colorbar;

                subplot(2,1,2);
                probeT = probe; probeT(badTimes,:) = [];
                [specT,~,~] = mtspecgramc(probeT(:,channels(channels>=15)),[5 2],params);
                meanNew = mean(10.*log10(abs(specT)),3,'omitnan'); % Mean across channels 15-32
                imagesc(timeValsSpec,freqValsSpec,meanNew'); colormap jet; shading interp; set(gca,'YDir','normal');
                xlabel('Time (s)'); ylabel('Frequency (Hz)'); clim([-30 30]); title('Spectrogram after removal of bad time segments'); colorbar; 

                sgtitle(strrep([expDate ' datafile ' num2str(fileNum) ' Probe: ' probeName ],'_','\_'));
                f = gcf; exportgraphics(f,[saveFolder '\AvgSpectrograms\AvgSpecgram_Before_After_' num2str(iFile) '_Probe' probeName '.png'],'Resolution',300);
                close gcf;
            end
           
        end

        allBadTimes{fileNum,iDate} = single(unique([badTimesA,badTimesB])); % Common bad time segments in both probes

        % Get intraprobe, interprobe correlograms and intraprobe verticals
        % for different frequencies

        chA = chInCortexProbeA{iDate}(iFile)+5; % 0.5mm from the transition point listed in experiment notes
        if ~strcmp(expDate,'08_08_2022')
            chB = chInCortexProbeB{iDate}(iFile)+5;
        else
            chB = 1;
        end

        for iBand = 1:3
            clear pA pB figTitle b a
            switch iBand
                case 1
                    pA = probe1;
                    pB = probe2;
                    figTitle = 'WB';

                case 2
                    pA = single(filtfilt(bG,aG,double(probe1)));
                    pB = single(filtfilt(bG,aG,double(probe2)));
                    figTitle = 'Gamma band';

                case 3
                    pA = single(filtfilt(bA,aA,double(probe1)));
                    pB = single(filtfilt(bA,aA,double(probe1)));
                    figTitle = 'Alpha band';
            end

            probeBList      = [1:13 15:21 23:32];
            if size(pA,2) == 33; pA(:,1) = []; end
            if size(pB,2) == 33; pB(:,1) = []; end
            if isempty(pA) || isempty(pB); continue; end

            % Remove bad channels
            if ~isempty(badElecA{fileNum,iDate} ); pA(:,badElecA{fileNum,iDate}) = []; end
            if ~isempty(badElecB{fileNum,iDate});  pB(:,badElecB{fileNum,iDate}) = [];  end


            % Remove the bad time segments
            if ~isempty(allBadTimes{fileNum,iDate})
                pA(allBadTimes{fileNum,iDate},:) = [];
                pB(allBadTimes{fileNum,iDate},:) = [];
            end

            clear cortexChA chOutA chDeepA cortexChB chOutB chDeepB
            cortexChA = pA(:,chA:chA+2);               % Correlating with channels that we know is inside the cortex
            chOutA    = pA(:,chOutCortex);             % Correlating with channels that we know is out of cortex
            chDeepA   = pA(:,size(pA,2)-2:size(pA,2)); % Correlating with the deepest channels

            if ~strcmp(expDate,'08_08_2022') % Condition where there is only one channel for a site (applicable for Whiskey)
                cortexChB = pB(:,chB:chB+2);
                chOutB    = pB(:,chOutCortex);
                chDeepB   = pB(:,size(pB,2)-2:size(pB,2));
            else
                cortexChB = pB;
                chOutB    = zeros(size(pB));
                chDeepB   = pB;
            end

            % Get the intra probe marginals...
            marginalA = mean(corr(pA,'Rows','complete'),2);
            marginalB = mean(corr(pB,'Rows','complete'),2);


            % Plot and save the intra probe correlation heatmaps...
            if ~exist([saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'file') ||  saveFigureFlag
                figure;
                imagesc(imgaussfilt(corr(pA,'Rows','complete'),1)); colormap jet; xticks(1:32); yticks(1:32); clim([ 0 1]); colorbar;
                title(strrep([expDate ' Datafile ' num2str(iFile) ' Probe A - ' figTitle ],'_','\_'));
                f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeA_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                close gcf;
            end

            if ~exist([saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'file')|| saveFigureFlag
                figure;
                imagesc(imgaussfilt(corr(pB,'Rows','complete'),1)); colormap jet; xticks(1:30); yticks(1:30); clim([ 0 1]); colorbar;
                xticklabels(probeBList); yticklabels(probeBList);title(strrep([expDate ' Datafile ' num2str(iFile) ' Probe B - ' figTitle  ],'_','\_'));
                f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\ProbeB_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                close gcf;
            end

            if ~exist([saveFolder '\IntraProbeCorr\Marginal_File_' num2str(fileNum) '_' figTitle '.png'],'file') || saveFigureFlag
                figure; subplot(1,2,1);

                imagesc(imgaussfilt(marginalA,1)); colormap jet; yticks(1:32); clim([0 1]);colorbar;  %caxis([ 0 max(marginalA)]);
                title('Probe A');

                subplot(1,2,2);
                imagesc(imgaussfilt(marginalB,1)); colormap jet; yticks(1:30);yticklabels(probeBList); clim([0 1]);

                colorbar;
                title('Probe B');
                sgtitle(strrep([expDate ' Datafile ' num2str(iFile) ' ' figTitle ],'_','\_'));

                f = gcf; exportgraphics(f,[saveFolder '\IntraProbeCorr\Marginal_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                close gcf;
            end

            % Plot and save the intra probe verticals...
            if ~exist([saveFolder '\Transition\Transition_' num2str(fileNum) '_Probe A_' figTitle '.png'],'file')  ||  saveFigureFlag
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
                                if ~strcmp(expDate,'08_08_2022')
                                    in = chB:chB+2;
                                else
                                    in = chB;
                                end
                                probeTitle = 'Probe B';
                        end

                        figure;
                        subplot(1,3,1); imagesc(movmean(mean(corrIn,2,'omitnan'),5));  clim([ 0 1]); colormap jet;
                        ylabel('Mean correlation'); title(['Ch: ' num2str(in)]);  xticks([]); yticks(1:size(corrIn,1)); colorbar;
                        if iProbe == 2; yticks(1:30); yticklabels(probeBList);end

                        subplot(1,3,2); imagesc(movmean(mean(corrOut,2,'omitnan'),5)); clim([ 0 1]); colormap jet;
                        ylabel('Mean correlation');title( ['Ch: ' num2str(chOutCortex)]);xticks([]); yticks(1:size(corrOut,1)); colorbar;
                        if iProbe == 2; yticks(1:30); yticklabels(probeBList);end

                        subplot(1,3,3); imagesc(movmean(mean(corrDeep,2,'omitnan'),5));    clim([0 1]); colormap jet;
                        ylabel('Mean correlation'); title(['Ch: ' num2str(chDeep)]);xticks([]); yticks(1:size(corrDeep,1)); colorbar;
                        if iProbe == 2; yticks(1:30); yticklabels(probeBList);end

                        sgtitle(strrep([expDate ' Datafile ' num2str(fileNum) ' ' probeTitle ' - ' figTitle  ],'_','\_'));

                        f = gcf; exportgraphics(f,[saveFolder '\Transition\Transition_' num2str(fileNum) '_' probeTitle '_' figTitle '.png'],'Resolution',300);
                        close gcf;
                    end
                end
            end

            % Plot and save the pairwise correlation heatmaps
            if ~exist([saveFolder '\InterProbeCorr\AllCh_File_' num2str(fileNum) '_' figTitle '.png'],'file') ||  saveFigureFlag
                figure; imagesc(imgaussfilt(corr(pA,pB,'rows','complete'),1)); colormap jet; clim([ 0 1]);
                title(strrep([expDate ' Datafile ' num2str(fileNum) ': Pairwise correlation -  ' figTitle],'_','\_')); xlabel('Probe B'); ylabel('Probe A'); colorbar;
                yticks(1:32); xticks(1:30); xticklabels(probeBList);
                f = gcf; exportgraphics(f,[saveFolder '\InterProbeCorr\AllCh_File_' num2str(fileNum) '_' figTitle '.png'],'Resolution',300);
                close gcf;
            end
        end
    end
end
end

