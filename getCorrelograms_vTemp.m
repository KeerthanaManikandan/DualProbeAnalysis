function [probe1,probe2,badElecA,badElecB] = getCorrelograms_vTemp(monkeyName, allDates, datFileNumAll,probe1,probe2,probeLabelA,probeLabelB,saveFigureFlag)

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
gammaBand       = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand       = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
badChProbeB     = [14 22];
probeBList      = [1:13 15:21 23:32];


for iDate =1:size(allDates,1)
    clc; clear expdate datFileNum datFileName
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results'];

    if ~exist([saveFolder '\AvgSpectrograms'],'dir'); [~,~] = mkdir([saveFolder '\AvgSpectrograms']); end
    if ~exist([saveFolder '\IntraProbeCorr'],'dir');  [~,~] = mkdir([saveFolder '\IntraProbeCorr']);  end
    if ~exist([saveFolder '\InterProbeCorr'],'dir');  [~,~] = mkdir([saveFolder '\InterProbeCorr']);  end
    if ~exist([saveFolder '\Transition'],'dir');      [~,~] = mkdir([saveFolder '\Transition']);      end

    for iFile = 1:length(datFileNum)
        if strcmp(expDate,'09_19_2022') && iFile == 4;  continue; end
        if strcmp(expDate,'02_07_2023') && iFile == 2; continue;  end
        if strcmp(expDate,'04_11_2023') && iFile == 10; continue; end

        fileNum = datFileNum(iFile);

        disp(['Preprocessing: ' num2str(expDate) ' File: ' num2str(fileNum)]);
        if isempty(probe1{fileNum,iDate}) || isempty(probe2{fileNum,iDate})
            badElecA{fileNum,iDate} = [];
            badElecB{fileNum,iDate} = [];
            continue;
        end
        probeA     = probe1{fileNum,iDate};
        channelsA  = 1:size(probeA,2);
        if ~strcmp(expDate,'08_08_2022') && ~isempty(probeLabelA{iDate,1})
            probeALabel = probeLabelA{iDate,1}(fileNum,:);
        else
            probeALabel = [];
        end
        probeAName = 'A';

        [spec,timeValsSpec,~] = mtspecgramc(probeA(:,channelsA),[5 2],params);
        powTimeBin = squeeze(sum(10.*log10(abs(spec)),2));
        if isempty(timeValsSpec); continue; end

        % Remove bad channel
        badElecA{fileNum,iDate} = [];
        if ~isempty(probeALabel) % Bad channel from previous
            if strcmp(probeALabel,'BD29')
                channelsA(ismember(channelsA,8)) = [];
                badElecA{fileNum,iDate} = [ badElecA{fileNum,iDate} 8];
            elseif strcmp(probeALabel,'CDE1')
                channelsA(ismember(channelsA,9)) = [];
                badElecA{fileNum,iDate} = [ badElecA{fileNum,iDate} 9];
            elseif strcmp(probeALabel,'D553')
                channelsA(ismember(channelsA,5)) = [];
                badElecA{fileNum,iDate} = [badElecA{fileNum,iDate} 5];
            end
        end
        badElecThresh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
        badElecA{fileNum,iDate} = [  badElecA{fileNum,iDate} channelsA(sum((powTimeBin>badElecThresh),1) >= floor(0.75*size(powTimeBin,1)))];
        if ~isempty( badElecA{fileNum,iDate}); channelsA(ismember(channelsA, badElecA{fileNum,iDate})) = []; end

        probeB     = probe2{fileNum,iDate};
        probeBName = 'B';
        if ~strcmp(expDate,'08_08_2022')
            channelsB = [1:13 15:21 23:size(probeB,2)];
            if ~isempty(probeLabelB{iDate,1})
                probeBLabel = probeLabelB{iDate,1}(fileNum,:);
            else
                probeBLabel = [];
            end
        else
            badTimesB = [];
            probeBLabel = [];
            badElecB{fileNum,iDate} = [];
            continue;
        end

        [spec,timeValsSpec,~] = mtspecgramc(probeB(:,channelsB),[5 2],params);
        powTimeBin = squeeze(sum(10.*log10(abs(spec)),2));
        if isempty(timeValsSpec); continue; end
        % Remove bad channel
        badElecB{fileNum,iDate} = [];
        if ~isempty(probeALabel) % Bad channel from previous
            if strcmp(probeALabel,'BD29')
                channelsA(ismember(channelsA,8)) = [];
                badElecB{fileNum,iDate} = [ badElecB{fileNum,iDate} 8];
            elseif strcmp(probeALabel,'CDE1')
                channelsA(ismember(channelsA,9)) = [];
                badElecB{fileNum,iDate} = [  badElecB{fileNum,iDate} 9];
            elseif strcmp(probeALabel,'D553')
                channelsA(ismember(channelsA,5)) = [];
              badElecB{fileNum,iDate} = [badElecB{fileNum,iDate} 5];
            end
        end
        badElecThresh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
        badElecB{fileNum,iDate} = [badElecB{fileNum,iDate} channelsB(sum((powTimeBin>badElecThresh),1) >= floor(0.75*size(powTimeBin,1)))];
        if ~isempty( badElecB{fileNum,iDate}); channelsB(ismember(channelsB, badElecB{fileNum,iDate})) = []; end


        % Finding bad time segments
        countVal = 1;  badTimeInd = [];
        while(countVal<=10)
            clear spec
            badTimes = [];
            if countVal<=5
                [spec,timeValsSpec,~] = mtspecgramc(probeA(:,channelsA(channelsA>=15)),[5 2],params);
            else
                [spec,timeValsSpec,~] = mtspecgramc(probeB(:,channelsB(channelsB>=15)),[5 2],params);
            end
            meanS = mean(10.*log10(abs(spec)),3,'omitnan'); % Mean across channels 15-32
            powMeanS = squeeze(sum(meanS,2));

            badTimeThresh   = (median(powMeanS,1)+4.5*mad(powMeanS,1,1));
            badTimeIndOld = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3;

            if ~isempty(badTimeIndOld)
                badTimeInd =[(badTimeIndOld-1e3)'  (badTimeIndOld+1e3)']; % Taking one second before and one second after bad time segments

                for iL = 1:size(badTimeInd,1)
                    badTimes = [ badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
                end
                probeA(badTimes,:) = [];
                probeB(badTimes,:) = [];
            end
            countVal = countVal+1;

        end

        probe1{fileNum,iDate} = probeA;
        probe2{fileNum,iDate} = probeB;
    end
end