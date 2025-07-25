%% Original data
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

for iDate = 1:size(allDates,1)
    clc; clear expdate datFileNum datFileName
    expDate = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    for iFile = 1:length(datFileNum)
        fileNum = datFileNum(iFile);
        disp(['Processing: ' num2str(expDate) ' File: ' num2str(fileNum)]); 
        for iProbe = 1:2
            clear probe channls spec timeVals powTimeBin badElecThresh badElecVal badTimeThresh badTimes
            if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate})
                badElecA{fileNum,iDate} = [];
                badElecB{fileNum,iDate} = [];
                continue; 
            end
            switch iProbe
                case 1
                    probe    = probe1{fileNum,iDate};
                    channels = 1:size(probe,2);
                case 2
                    probe = probe2{fileNum,iDate};
                    if ~strcmp(expDate,'08_08_2022')
                    channels = [1:13 15:21 23:size(probe,2)];
                    else
                        badTimesB = [];
                        badElecB{fileNum,iDate} = [];
                        continue;
                    end
            end
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probe(:,channels),[5 2],params);
            powTimeBin = squeeze(sum(10.*log10(abs(spec)),2));  
            if isempty(timeValsSpec); continue; end

            % Remove bad electrode
            badElecThresh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
            badElecVal = channels(sum((powTimeBin>badElecThresh),1) >= floor(0.75*size(powTimeBin,1)));
           
            % Remove bad time points
            if ~isempty(badElecVal); channels(ismember(channels,badElecVal)) = []; end 
            clear spec
            [spec,~,~] = mtspecgramc(probe(:,channels),[5 2],params);
            meanS = mean(10.*log10(abs(spec)),3,'omitnan');
            powMeanS = squeeze(sum(meanS,2));

            badTimeThresh   = (median(powMeanS,1)+5*mad(powMeanS,1,1));
            badTimeInd = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3;
             badTimes = [];
            if ~isempty(badTimeInd)
                l = find(isoutlier([0 diff(badTimeInd)]));
                for iL = 1:length(l)
                    if iL == length(l)
                        badTimes = [ badTimes badTimeInd(l(iL)): badTimeInd(end)];
                    else

                        badTimes = [ badTimes badTimeInd(l(iL)): badTimeInd(l(iL+1)-1)];
                    end
                end
            end

            switch iProbe
                case 1
                    badTimesA               = badTimes;
                    badElecA{fileNum,iDate} = badElecVal;
                case 2
                    badTimesB               = badTimes;
                    badElecB{fileNum,iDate} = badElecVal;
            end
        end

        allBadTimes{fileNum,iDate} = unique([badTimesA,badTimesB]);

        % Find channels outside of cortex using PCA
        cVals = [1 0 0 ; 0 0 1];
        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate})
            transitionVals{fileNum,iDate} = [];
            clusterIdx{fileNum,iDate}=[];
            continue;
        end
        for iProbe = 1:2
            clear powTimeBin probe channels figtitle 
          
            switch iProbe
                case 1
                    probe = probe1{fileNum,iDate};
                    if size(probe,2) == 33; probe(:,1) = []; end
                    channels = 1:size(probe,2);
                    badElec = badElecA{fileNum,iDate};
                    figTitle = 'Probe A';
                case 2
                    probe = probe2{fileNum,iDate};
                    if size(probe,2) == 33; probe(:,1) = []; end
                    badElec = badElecB{fileNum,iDate};
                    if ~strcmp(expDate,'08_08_2022')
                        channels = [1:13 15:21 23:size(probe,2)];
                    else
                        transitionVals{fileNum,iDate}{iProbe} = 1;
                        clusterIdx{fileNum,iDate}{iProbe} =[];
                        continue;
                    end
                    figTitle = 'Probe B';
            end
           
            % Remove bad time segments and bad electrodes
            clear spec   
            try probe(allBadTimes{fileNum,iDate},:);
                probe(allBadTimes{fileNum,iDate},:) = [];
            catch
                continue;
            end
            if ~isempty(badElec); channels(ismember(channels,badElec)) = []; end  % Remove bad electrodes

           if isempty(probe)
                transitionVals{fileNum,iDate}{iProbe} = 0;
                clusterIdx{fileNum,iDate}{iProbe} = [];
                continue;
           end

            [spec,~,freqVals] = mtspecgramc(probe(:,channels),[5 2],params);
%             pp = squeeze(sum(10.*log10(abs(spec)),2));
            
            if isempty(spec)
                transitionVals{fileNum,iDate}{iProbe} = 0;
                clusterIdx{fileNum,iDate}{iProbe} = [];
                continue;
            end
            freqIndex = find(freqVals>=20 & freqVals<=55);
            powTimeBin = squeeze(sum(10.*log10(abs(spec(:,freqIndex,:))),2));

            % Perform PCA 
            clear mu N S Um lambda PC_scores idx
            mu = mean(powTimeBin,2);
            N = size(powTimeBin,2);
            S = (bsxfun(@minus,powTimeBin,mu)*bsxfun(@minus,powTimeBin,mu)')./N; % Co-variance matrix
            [Um,lambda] = eig(S);
            PC_scores = Um(:,1:2)'*bsxfun(@minus,powTimeBin,mu);

            % Perform k-means to separate the two clusters
            idx = kmeans(PC_scores',2);
            idxDiff = find([0; abs(diff(idx))]);
            idxDiff(idxDiff>=20) = [];
%             [~,t] = min(16-idxDiff);
            transitionVals{fileNum,iDate}{iProbe} = median(idxDiff);
            clusterIdx{fileNum,iDate}{iProbe} = idx;
%             figure;scatter(PC_scores(1,:)',PC_scores(2,:)',30,cVals(idx,:),'filled');
%             grid on; xlabel('PC 1 score'); ylabel('PC 2 score'); box on;
%             title(figTitle);
        end
    end
end


% [coeff,score,latent,tsq,explained,mu] = pca(ssANew');

%%
rowIdx = 1; 
bandLabels = {'Alpha band'; 'Beta band'; 'Gamma band';'Wideband'}; 
tic;
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
            corrVals(rowIdx,iBand) = NaN;
            rowIdx = rowIdx+1; continue; 
        end 
        try probeA(allBadTimes{fileNum,iDate},:);
            probeA(allBadTimes{fileNum,iDate},:) = [];
        catch
            corrVals(rowIdx,iBand) = NaN;
            rowIdx = rowIdx+1;
            continue;
        end

        try probeB(allBadTimes{fileNum,iDate},:);
            probeB(allBadTimes{fileNum,iDate},:) = [];
        catch
            corrVals(rowIdx,iBand) = NaN;
            rowIdx = rowIdx+1;
            continue;
        end
%         probeA(allBadTimes{fileNum,iDate},:) = [];
%         probeB(allBadTimes{fileNum,iDate},:) = [];

        if ~isempty(badElecA); probeA(:,badElecA{fileNum,iDate}) = []; end 
        if ~isempty(badElecB); probeA(:,badElecB{fileNum,iDate}) = []; end


        if isempty(probeA) || isempty(probeB)
            corrVals(rowIdx,iBand) = NaN;
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
                corrVals(rowIdx,iBand) = NaN;
                continue;
            end

            if isempty(probe1{fileNum,iDate}) && isempty(probe2{fileNum,iDate})
                corrVals(rowIdx,iBand) = NaN;
                continue;
            end

            clear xA yA chA chB
            if isempty(transitionVals{fileNum,iDate}{1}) || isempty(transitionVals{fileNum,iDate}{2})
                corrVals(rowIdx,iBand) = NaN;
                continue;
            end

            chA(1) = transitionVals{fileNum,iDate}{1};
            if chA(1)+19 >size(probeA,2)
                chA(2) = size(probeA,2);
            else
                chA(2) = chA(1)+19;
            end 
            if chA(1) == 0 || isnan(chA(1))
                corrVals(rowIdx,iBand) = NaN;
                continue;
            end
            
            chB(1) = transitionVals{fileNum,iDate}{2};

            if chB(1) == 0 || isnan(chB(1))
                 corrVals(rowIdx,iBand) = NaN;
                continue;
            end

            if chB(1)+19 >size(probeB,2)
                chB(2) = size(probeB,2);
            else
                chB(2) = chB(1)+19;
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
            corrVals(rowIdx,iBand) = mean(corr(xA,yA),'all','omitnan'); % Taking the mean correlations
        end 
         rowIdx = rowIdx+1;
    end
   
end
nanVals = find(isnan(corrVals(:,4)));
corrVals(nanVals,:) = [];
%% Run this once only
connValsAll(nanVals) = []; % 
%% Plotting
for iBand = 1:4
  
    subplot(2,2,iBand);
    scatter(connValsAll,squeeze(corrVals(:,iBand)),20,'filled');
    xlabel('Functional connectivity'); xlim([-0.3 1]);
    ylabel('Pairwise correlations'); ylim([-1 1]);
    hold on; box on;title(bandLabels{iBand}); grid on;

    xFit = linspace(min(connValsAll),max(connValsAll),1000);
    [f,gof]  = fit(connValsAll,squeeze(corrVals(:,iBand)),'poly1');

    c        = coeffvalues(f);
    r2       = gof.rsquare;
    fitLine   = c(1)*(xFit) +c(2);
    plot(xFit,fitLine,'Color','k','LineStyle','--','LineWidth',2);
    text(0.6,-0.8,['R^2: ' num2str(r2)]);
end 
