function [removeTimes] = getCorrelograms(monkeyName, allDates, datFileNumAll,probe1,probe2,chInCortexProbeA,chInCortexProbeB,saveFigureFlag)

if ~exist('saveFigureFlag','var') || isempty(saveFigureFlag); saveFigureFlag = 0; end 
hemisphere  = 'Left';
removeTimes = cell(size(probe1)); 
badChProbeB = [14 22];
fs          = 1e3; % Sampling frequency
chOutCortex = 1:3;
chDeep      = 30:32;
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand   = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
probeBList  = [1:13 15:21 23:32];

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

end