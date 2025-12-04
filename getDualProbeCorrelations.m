function [allVars] = getDualProbeCorrelations(monkeyName, hemisphere, allDates, datFileNumAll,allProbeData,allBadTimes,...
    badElecA,badElecB,estChInCortexA,estChInCortexB,connValsAll,distSitesAll,goodRuns)
clear meanPSDEEG medPSDEEG  specAMeanAll specBMeanAll specAMedAll specBMedAll meanCorrInfraSlow medCorrInfraSlow
fs = 1e3;
[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);
                            
thetaBand = [6 8];    [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters
alphaBand = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass');
betaBand  = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass');

bandLabels = {'Theta', 'Alpha', 'Beta', 'Gamma','Spiking'};
freqCombs = nchoosek(1:size(bandLabels,2),2);

% Check if certain variables are stored....
if exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') 
    varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    varIdx  = ismember('infraIntraAAllR',varInfo);
else 
    varIdx = 1;
end

% Calculate pairwise correlations between and within probes
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'file') || ~varIdx
    tic; % re-run for Whiskey... 
    for iDate = 1: size(allDates,1)
        clear expDate datFileNum saveFolder
        expDate    = allDates(iDate,:);
        datFileNum = datFileNumAll{iDate,1};

        for iRun = 1:length(datFileNum)
            fileNum = datFileNum(iRun);
            clear probeA probeB chA chB

            clc; disp(['Processing data for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(fileNum)]);

            % Get the ephys data
            clear probeA probeB rawA rawB
            probeA = allProbeData{fileNum,iDate}.probe1Ch;
            probeB = allProbeData{fileNum,iDate}.probe2Ch;

            rawA   = allProbeData{fileNum,iDate}.raw1Ch;
            rawB   = allProbeData{fileNum,iDate}.raw2Ch;

            % Remove bad channels
            probeA(:,badElecA{fileNum,iDate}) = [];
            probeB(:,badElecB{fileNum,iDate}) = [];

            rawA(:,badElecA{fileNum,iDate}) = [];
            rawB(:,badElecB{fileNum,iDate}) = [];

            % Remove bad times
            probeA(allBadTimes{fileNum,iDate},:) = [];
            probeB(allBadTimes{fileNum,iDate},:) = [];

            rawA(allBadTimes{fileNum,iDate},:) = [];
            rawB(allBadTimes{fileNum,iDate},:) = [];

            chA = estChInCortexA{iDate}(iRun,:);
            chB = estChInCortexB{iDate}(iRun,:);

            if chA(1) == 0 || chB(1) == 0 || isempty(probeA) || isempty(probeB)                
                meanPairCorr(iDate,iRun,1:5)   = NaN; maxPairCorr(iDate,iRun,1:5)    = NaN;
                medPairCorr(iDate,iRun,1:5)    = NaN; meanIntraCorrA(iDate,iRun,1:5) = NaN;
                meanIntraCorrB(iDate,iRun,1:5) = NaN; medIntraCorrA(iDate,iRun,1:5)  = NaN;
                medIntraCorrB(iDate,iRun,1:5)  = NaN;

                medPairCorrSuper(iDate,iRun,1:5)      = NaN;
                medPairCorrMid(iDate,iRun,1:5)        = NaN;
                medPairCorrDeep(iDate,iRun,1:5)       = NaN;

                envelopePairCorrSuper(iDate,iRun,1:5) = NaN;
                envelopePairCorrMid(iDate,iRun,1:5)   = NaN;
                envelopePairCorrDeep(iDate,iRun,1:5)  = NaN;

                infraPairCorrSuper(iDate,iRun,1:5)    = NaN;
                infraPairCorrMid(iDate,iRun,1:5)      = NaN;
                infraPairCorrDeep(iDate,iRun,1:5)     = NaN;

                medPairCorrFreqTime(1:size(freqCombs,2),iDate,iRun)= NaN;   % Cross frequency pairs
                medPairCorrFreqPower(1:size(freqCombs,2),iDate,iRun)= NaN;
                medPairCorrFreqInfra(1:size(freqCombs,2),iDate,iRun)= NaN;

                continue;
            end

            % Get mean, median and maximum pairwise correlations for different
            % frequency bands...
            clear timeSeriesCorr sizeSpecA sizeSpecB specA_R specB_R freqSpecCorr
            for iBand = 1:5
                clear xA xB xARef xBRef specA specB specEEG timeValsSpec freqValsSpec eegGood envelopeABandLimited envelopeBBandLimited infraSlowA infraSlowB

                if (strcmp(expDate,'10_10_2022') && (fileNum == 1 || fileNum == 6 ))||...
                        (strcmp(expDate,'02_07_2023') && (fileNum == 2 ))... %||(strcmp(expDate,'11_01_2021') && (fileNum == 11 ))...
                        ||(strcmp(expDate,'09_19_2022') && (fileNum == 4)) || (isempty(probeA) && isempty(probeB))...
                        || (strcmp(expDate,'04_11_2023') && (fileNum == 10))% Datafile 4 from 09/19/2022 - why -  (strcmp(expDate,'09_19_2022') && fileNum == 4) ||

                    meanPairCorr(iDate,iRun,iBand)   = NaN; maxPairCorr(iDate,iRun,iBand)    = NaN; medPairCorr(iDate,iRun,iBand)  = NaN;
                    meanIntraCorrA(iDate,iRun,iBand) = NaN; meanIntraCorrB(iDate,iRun,iBand) = NaN;
                    medIntraCorrA(iDate,iRun,iBand)  = NaN; medIntraCorrB(iDate,iRun,iBand)  = NaN;

                    medPairCorrSuper(iDate,iRun,iBand)  = NaN; medPairCorrMid(iDate,iRun,iBand)  = NaN; medPairCorrDeep(iDate,iRun,iBand)  = NaN;
                    envelopePairCorrSuper(iDate,iRun,iBand)  = NaN; envelopePairCorrMid(iDate,iRun,iBand)  = NaN; envelopePairCorrDeep(iDate,iRun,iBand)  = NaN;
                    infraPairCorrSuper(iDate,iRun,iBand)  = NaN; infraPairCorrMid(iDate,iRun,iBand)  = NaN; infraPairCorrDeep(iDate,iRun,iBand)  = NaN;

                    continue;
                end

                if chA(1) == 0 || chB(1)== 0 % Check if any of the recordings are empty
                    meanPairCorr(iDate,iRun,iBand)   = NaN; maxPairCorr(iDate,iRun,iBand)    = NaN;
                    medPairCorr(iDate,iRun,iBand)    = NaN; meanIntraCorrA(iDate,iRun,iBand) = NaN;
                    meanIntraCorrB(iDate,iRun,iBand) = NaN; medIntraCorrA(iDate,iRun,iBand)  = NaN;
                    medIntraCorrB(iDate,iRun,iBand)  = NaN;
                    continue;
                end

                clear xA yA
                switch iBand % Get the correlation vs connectivity vs distance plots for different bands...
                    case 1 % Theta band
                        xA = filtfilt(bT,aT,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bT,aT,double(probeB(:,chB(1):chB(2))));

                    case 2 % Alpha band
                        xA = filtfilt(bA,aA,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bA,aA,double(probeB(:,chB(1):chB(2))));

                    case 3 % Beta band
                        xA = filtfilt(bB,aB,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bB,aB,double(probeB(:,chB(1):chB(2))));

                    case 4 % Gamma band
                        xA = filtfilt(bG,aG,double(probeA(:,chA(1):chA(2))));
                        xB = filtfilt(bG,aG,double(probeB(:,chB(1):chB(2))));

                    case 5 % Spiking
                        xA = rawA(:,chA(1):chA(2));
                        xB = rawB(:,chB(1):chB(2));
                end

                % Get the intra probe correlations for channels inside the
                % cortex...
                medIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(xA,'rows','complete'))),'all','omitnan');
                medIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(xB,'rows','complete'))),'all','omitnan');
               
                intraCorrA{iDate,iRun,iBand} = single(tril(corr(xA,'rows','complete')));
                intraCorrB{iDate,iRun,iBand} = single(tril(corr(xB,'rows','complete')));

                % Get pairwise correlations between the two probes...
                maxPairCorr(iDate,iRun,iBand)  = max(corr(xA,xB),[],'all','omitnan');
                meanPairCorr(iDate,iRun,iBand) = mean(corr(xA,xB),'all','omitnan');
                medPairCorr(iDate,iRun,iBand)  = median(corr(xA,xB),'all','omitnan');%

                % Get instantaneous power and correlate the powers
                envelopeABandLimited = envelope(abs(xA),5);
                envelopeBBandLimited = envelope(abs(xB),5);

                % Within probe correlations - Envelope
                envelopeIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(envelopeABandLimited,envelopeABandLimited))),'all','omitnan');
                envelopeIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(envelopeBBandLimited,envelopeBBandLimited))),'all','omitnan');

                envelopeIntraAAll{iDate,iRun,iBand} =  single(tril(corr(envelopeABandLimited,'rows','complete')));
                envelopeIntraBAll{iDate,iRun,iBand} =  single(tril(corr(envelopeBBandLimited,'rows','complete')));

                % Correlate instantaneous band power
                medCorrEnvelope(iDate,iRun,iBand)  =  median(corr(envelopeABandLimited,envelopeBBandLimited),'all','omitnan');

                enSizeA = size(envelopeABandLimited);
                enSizeB = size(envelopeBBandLimited);

                % Filter envelope from 0.01-0.1 Hz
                envelopeABandLimited = [envelopeABandLimited; envelopeABandLimited ;envelopeABandLimited];
                infraSlowA = filtfilt(sos,g,double(envelopeABandLimited));
                infraSlowA = single(infraSlowA(enSizeA(1)+1:(end-enSizeA(1)),:));

                envelopeBBandLimited = [envelopeBBandLimited; envelopeBBandLimited ;envelopeBBandLimited];
                infraSlowB = filtfilt(sos,g,double(envelopeBBandLimited));
                infraSlowB = single(infraSlowB(enSizeB(1)+1:(end-enSizeB(1)),:));

                % Within probe correlations - Infraslow
                infraIntraCorrA(iDate,iRun,iBand)  = median(nonzeros(tril(corr(infraSlowA,infraSlowA))),'all','omitnan');
                infraIntraCorrB(iDate,iRun,iBand)  = median(nonzeros(tril(corr(infraSlowB,infraSlowB))),'all','omitnan');

                infraIntraAAll{iDate,iRun,iBand} =  single(tril(corr(infraSlowA,'rows','complete')));
                infraIntraBAll{iDate,iRun,iBand} =  single(tril(corr(infraSlowB,'rows','complete')));

                % Correlate infraslow flucutuations in instantaneous band power
                medCorrInfraSlow(iDate,iRun,iBand)  =  median(corr(infraSlowA,infraSlowB),'all','omitnan');


                % Split the data into superficial, middle or deep for both channels
                chSplit = 6;
                % Electrode A
                if chA(2)-chA(1)== 0 % Single channel
                    timeASuper = xA;
                    timeAMid   = xA;
                    timeADeep  = xA;

                    envelopeASuper = envelopeABandLimited;
                    envelopeAMid   = envelopeABandLimited;
                    envelopeADeep  = envelopeABandLimited;

                    infraASuper = infraSlowA;
                    infraAMid   = infraSlowA;
                    infraADeep  = infraSlowA;
                else
                    timeASuper = xA(:,1:chSplit);
                    timeAMid   = xA(:,(chSplit+1: chSplit*2));
                    timeADeep  = xA(:,2*chSplit+1:end);

                    envelopeASuper = envelopeABandLimited(:,1:chSplit);
                    envelopeAMid   = envelopeABandLimited(:,(chSplit+1: chSplit*2));
                    envelopeADeep  = envelopeABandLimited(:,2*chSplit+1:end);

                    infraASuper = infraSlowA(:,1:chSplit);
                    infraAMid   = infraSlowA(:,(chSplit+1:chSplit*2));
                    infraADeep  = infraSlowA(:,2*chSplit+1:end);

                end

                % Electrode B
                if chB(2)-chB(1)== 0 % Single channel
                    timeBSuper = xB;
                    timeBMid   = xB;
                    timeBDeep  = xB;

                    envelopeBSuper = envelopeBBandLimited;
                    envelopeBMid   = envelopeBBandLimited;
                    envelopeBDeep  = envelopeBBandLimited;

                    infraBSuper = infraSlowB;
                    infraBMid   = infraSlowB;
                    infraBDeep  = infraSlowB;
                else
                    timeBSuper = xB(:,1:chSplit);%
                    timeBMid   = xB(:,(chSplit+1: chSplit*2));%
                    timeBDeep  = xB(:,2*chSplit+1:end);%

                    envelopeBSuper = envelopeBBandLimited(:,1:chSplit);%
                    envelopeBMid   = envelopeBBandLimited(:,(chSplit+1: chSplit*2));%
                    envelopeBDeep  = envelopeBBandLimited(:,2*chSplit+1:end);%

                    infraBSuper = infraSlowB(:,1:chSplit); %
                    infraBMid   = infraSlowB(:,(chSplit+1: chSplit*2));%
                    infraBDeep  = infraSlowB(:,2*chSplit+1:end);%
                end

                % Correlate within compartments
                medPairCorrSuper(iDate,iRun,iBand) = median(corr(timeASuper,timeBSuper,'rows','complete'),'all','omitnan');
                medPairCorrMid(iDate,iRun,iBand)   = median(corr(timeAMid,timeBMid,'rows','complete'),'all','omitnan');
                medPairCorrDeep(iDate,iRun,iBand)  = median(corr(timeADeep,timeBDeep,'rows','complete'),'all','omitnan');

                envelopePairCorrSuper(iDate,iRun,iBand) = median(corr(envelopeASuper,envelopeBSuper,'rows','complete'),'all','omitnan');
                envelopePairCorrMid(iDate,iRun,iBand)   = median(corr(envelopeAMid,envelopeBMid,'rows','complete'),'all','omitnan');
                envelopePairCorrDeep(iDate,iRun,iBand)  = median(corr(envelopeADeep,envelopeBDeep,'rows','complete'),'all','omitnan');

                infraPairCorrSuper(iDate,iRun,iBand) = median(corr(infraASuper,infraBSuper,'rows','complete'),'all','omitnan');
                infraPairCorrMid(iDate,iRun,iBand)   = median(corr(infraAMid,infraBMid,'rows','complete'),'all','omitnan');
                infraPairCorrDeep(iDate,iRun,iBand)  = median(corr(infraADeep,infraBDeep,'rows','complete'),'all','omitnan');

                % Correlate between compartments between probes
                medPairASuperBMid(iDate,iRun,iBand)  = median(corr(timeASuper,timeBMid,'rows','complete'),'all','omitnan');
                medPairASuperBDeep(iDate,iRun,iBand) = median(corr(timeASuper,timeBDeep,'rows','complete'),'all','omitnan');
                medPairAMidBSuper(iDate,iRun,iBand)  = median(corr(timeAMid,timeBSuper,'rows','complete'),'all','omitnan');
                medPairAMidBDeep(iDate,iRun,iBand)   = median(corr(timeAMid,timeBDeep,'rows','complete'),'all','omitnan');
                medPairADeepBSuper(iDate,iRun,iBand) = median(corr(timeADeep,timeBSuper,'rows','complete'),'all','omitnan');
                medPairADeepBMid(iDate,iRun,iBand)   = median(corr(timeADeep,timeBMid,'rows','complete'),'all','omitnan');

                envelopeASuperBMid(iDate,iRun,iBand)  = median(corr(envelopeASuper,envelopeBMid,'rows','complete'),'all','omitnan');
                envelopeASuperBDeep(iDate,iRun,iBand) = median(corr(envelopeASuper,envelopeBDeep,'rows','complete'),'all','omitnan');
                envelopeAMidBSuper(iDate,iRun,iBand)  = median(corr(envelopeAMid,envelopeBSuper,'rows','complete'),'all','omitnan');
                envelopeAMidBDeep(iDate,iRun,iBand)   = median(corr(envelopeAMid,envelopeBDeep,'rows','complete'),'all','omitnan');
                envelopeADeepBSuper(iDate,iRun,iBand) = median(corr(envelopeADeep,envelopeBSuper,'rows','complete'),'all','omitnan');
                envelopeADeepBMid(iDate,iRun,iBand)   = median(corr(envelopeADeep,envelopeBMid,'rows','complete'),'all','omitnan');

                infraASuperBMid(iDate,iRun,iBand)  = median(corr(infraASuper,infraBMid,'rows','complete'),'all','omitnan');
                infraASuperBDeep(iDate,iRun,iBand) = median(corr(infraASuper,infraBDeep,'rows','complete'),'all','omitnan');
                infraAMidBSuper(iDate,iRun,iBand)  = median(corr(infraAMid,infraBSuper,'rows','complete'),'all','omitnan');
                infraAMidBDeep(iDate,iRun,iBand)   = median(corr(infraAMid,infraBDeep,'rows','complete'),'all','omitnan');
                infraADeepBSuper(iDate,iRun,iBand) = median(corr(infraADeep,infraBSuper,'rows','complete'),'all','omitnan');
                infraADeepBMid(iDate,iRun,iBand)   = median(corr(infraADeep,infraBMid,'rows','complete'),'all','omitnan');

                % Get cross-frequency correlations
                for iComb = 1:size(freqCombs,1)
                    clear combVals pAFreq pBFreq envelopeABand envelopeBBand
                    combVals = freqCombs(iComb,:);

                    for iPair = 1:length(combVals)
                        switch combVals(iPair)
                            case 1
                                pAFreq(:,:,iPair) = single(filtfilt(bT,aT,double(probeA(:,chA(1):chA(2)))));
                                pBFreq(:,:,iPair) = single(filtfilt(bT,aT,double(probeB(:,chB(1):chB(2)))));
                            case 2
                                pAFreq(:,:,iPair) = single(filtfilt(bA,aA,double(probeA(:,chA(1):chA(2)))));
                                pBFreq(:,:,iPair) = single(filtfilt(bA,aA,double(probeB(:,chB(1):chB(2)))));
                            case 3
                                pAFreq(:,:,iPair) = single(filtfilt(bB,aB,double(probeA(:,chA(1):chA(2)))));
                                pBFreq(:,:,iPair) = single(filtfilt(bB,aB,double(probeB(:,chB(1):chB(2)))));
                            case 4
                                pAFreq(:,:,iPair) = single(filtfilt(bG,aG,double(probeA(:,chA(1):chA(2)))));
                                pBFreq(:,:,iPair) = single(filtfilt(bG,aG,double(probeB(:,chB(1):chB(2)))));
                            case 5
                                pAFreq(:,:,iPair) = rawA(:,chA(1):chA(2));
                                pBFreq(:,:,iPair) = rawB(:,chB(1):chB(2));
                        end

                        % Get instantaneous power and correlate the powers
                        envelopeABand(:,:,iPair) = envelope(abs(pAFreq(:,:,iPair)),5);
                        envelopeBBand(:,:,iPair) = envelope(abs(pBFreq(:,:,iPair)),5);
                    end

                    medPairCorrFreqTime(iComb,iDate,iRun) = median([corr(squeeze(pAFreq(:,:,1)),squeeze(pBFreq(:,:,2))) corr(squeeze(pAFreq(:,:,2)),squeeze(pBFreq(:,:,1)))],'all','omitnan');
                    medPairCorrFreqTimeIntraA(iComb,iDate,iRun) = median(corr(squeeze(pAFreq(:,:,1)),squeeze(pAFreq(:,:,2))),'all','omitnan');
                    medPairCorrFreqTimeIntraB(iComb,iDate,iRun) = median(corr(squeeze(pBFreq(:,:,1)),squeeze(pBFreq(:,:,2))),'all','omitnan');

                    % Correlate instantaneous band power
                    medPairCorrFreqPower(iComb,iDate,iRun) = median([corr(squeeze(envelopeABand(:,:,1)),squeeze(envelopeBBand(:,:,2))) corr(squeeze(envelopeABand(:,:,2)),squeeze(envelopeBBand(:,:,1)))],'all','omitnan');
                    medPairCorrFreqPowerIntraA(iComb,iDate,iRun) = median(corr(squeeze(envelopeABand(:,:,1)),squeeze(envelopeABand(:,:,2))) ,'all','omitnan');
                    medPairCorrFreqPowerIntraB(iComb,iDate,iRun) = median(corr(squeeze(envelopeBBand(:,:,1)),squeeze(envelopeBBand(:,:,2))) ,'all','omitnan');

                    enSizeA = size(envelopeABand);
                    enSizeB = size(envelopeBBand);

                    % Filter envelope from 0.01-0.1 Hz
                    envelopeABand = [envelopeABand; envelopeABand ;envelopeABand];
                    infraSlowA = filtfilt(sos,g,double(envelopeABand));
                    infraSlowA = single(infraSlowA(enSizeA(1)+1:(end-enSizeA(1)),:,:));

                    envelopeBBand = [envelopeBBand; envelopeBBand ;envelopeBBand];
                    infraSlowB = filtfilt(sos,g,double(envelopeBBand));
                    infraSlowB = single(infraSlowB(enSizeB(1)+1:(end-enSizeB(1)),:,:));

                    % Correlate infraslow flucutuations in instantaneous band power
                    medPairCorrFreqInfra(iComb,iDate,iRun) =  median([corr(squeeze(infraSlowA(:,:,1)),squeeze(infraSlowB(:,:,2))) corr(squeeze(infraSlowA(:,:,2)),squeeze(infraSlowB(:,:,1)))],'all','omitnan');
                    medPairCorrFreqInfraIntraA(iComb,iDate,iRun) =  median(corr(squeeze(infraSlowA(:,:,1)),squeeze(infraSlowA(:,:,2))),'all','omitnan');
                    medPairCorrFreqInfraIntraB(iComb,iDate,iRun) =  median(corr(squeeze(infraSlowB(:,:,1)),squeeze(infraSlowB(:,:,2))),'all','omitnan');
                end

            end
        end
    end
    

    % Reshape variables for storage
    disp('Reshaping all variables to store as mat files...')
    matSize   = size(connValsAll);
    connValsR = reshape(connValsAll',[matSize(1)*matSize(2) 1]);
    distValsR = reshape(distSitesAll',[matSize(1)*matSize(2) 1]);

    intraCorrAR    = reshape(intraCorrA,[matSize(1)*matSize(2) size(intraCorrA,3)]);
    intraCorrBR    = reshape(intraCorrB,[matSize(1)*matSize(2) size(intraCorrB,3)]);

    envelopeIntraAAllR = reshape(envelopeIntraAAll,[matSize(1)*matSize(2) size(intraCorrA,3)]);
    envelopeIntraBAllR = reshape(envelopeIntraBAll,[matSize(1)*matSize(2) size(intraCorrA,3)]);

    infraIntraAAllR  = reshape(infraIntraAAll,[matSize(1)*matSize(2) size(intraCorrA,3)]); 
    infraIntraBAllR  = reshape(infraIntraAAll,[matSize(1)*matSize(2) size(intraCorrA,3)]);

    medPairCorrR      = reshape(medPairCorr,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medCorrEnvelopeR  = reshape(medCorrEnvelope,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    medCorrInfraSlowR = reshape(medCorrInfraSlow,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    medIntraCorrAR    = reshape(medIntraCorrA,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    medIntraCorrBR    = reshape(medIntraCorrB,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    envelopeIntraCorrAR = reshape(envelopeIntraCorrA,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    envelopeIntraCorrBR = reshape(envelopeIntraCorrB,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    infraIntraCorrAR = reshape(infraIntraCorrA,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    infraIntraCorrBR = reshape(infraIntraCorrB,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    medPairCorrSuperR = reshape(medPairCorrSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairCorrMidR   = reshape(medPairCorrMid,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    medPairCorrDeepR  = reshape(medPairCorrDeep,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    envelopePairCorrSuperR = reshape(envelopePairCorrSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopePairCorrMidR   = reshape(envelopePairCorrMid,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    envelopePairCorrDeepR  = reshape(envelopePairCorrDeep,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    infraPairCorrSuperR = reshape(infraPairCorrSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraPairCorrMidR   = reshape(infraPairCorrMid,[matSize(1)*matSize(2) size(meanPairCorr,3)]);
    infraPairCorrDeepR  = reshape(infraPairCorrDeep,[matSize(1)*matSize(2) size(meanPairCorr,3)]);

    medPairASuperBMidR  = reshape(medPairASuperBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairASuperBDeepR = reshape(medPairASuperBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairAMidBSuperR  = reshape(medPairAMidBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairAMidBDeepR   = reshape(medPairAMidBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairADeepBSuperR = reshape(medPairADeepBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    medPairADeepBMidR   = reshape(medPairADeepBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);

    envelopeASuperBMidR  = reshape(envelopeASuperBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeASuperBDeepR = reshape(envelopeASuperBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeAMidBSuperR  = reshape(envelopeAMidBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeAMidBDeepR   = reshape(envelopeAMidBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeADeepBSuperR = reshape(envelopeADeepBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    envelopeADeepBMidR   = reshape(envelopeADeepBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);

    infraASuperBMidR  = reshape(infraASuperBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraASuperBDeepR = reshape(infraASuperBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraAMidBSuperR  = reshape(infraAMidBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraAMidBDeepR   = reshape(infraAMidBDeep,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraADeepBSuperR = reshape(infraADeepBSuper,[matSize(1)*matSize(2) size(medPairCorr,3)]);
    infraADeepBMidR   = reshape(infraADeepBMid,[matSize(1)*matSize(2) size(medPairCorr,3)]);

    medPairCorrFreqTimeR        = reshape(medPairCorrFreqTime,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqTimeIntraAR  = reshape(medPairCorrFreqTimeIntraA,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqTimeIntraBR  = reshape(medPairCorrFreqTimeIntraB,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqPowerR       = reshape(medPairCorrFreqPower,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqPowerIntraAR = reshape(medPairCorrFreqPowerIntraA,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqPowerIntraBR = reshape(medPairCorrFreqPowerIntraB,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqInfraR       = reshape(medPairCorrFreqInfra,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqInfraIntraAR = reshape(medPairCorrFreqInfraIntraA,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);
    medPairCorrFreqInfraIntraBR = reshape(medPairCorrFreqInfraIntraB,[size(medPairCorrFreqInfra,1) size(medPairCorrFreqInfra,2)*size(medPairCorrFreqInfra,3)]);


    nanVals       = (isnan(connValsR) | isnan(medPairCorrR(:,1)));
    lowCorrVals   = (medIntraCorrAR(:,4)<=0.2 | medIntraCorrBR(:,4)<=0.2);
    removeDataIdx = nanVals | lowCorrVals;

    connValsR(removeDataIdx)             = [];
    distValsR(removeDataIdx)             = [];
    medPairCorrR(removeDataIdx,:)        = [];
    medCorrEnvelopeR(removeDataIdx,:)    = [];
    medCorrInfraSlowR(removeDataIdx,:)   = [];
    medIntraCorrAR(removeDataIdx,:)      = [];
    medIntraCorrBR(removeDataIdx,:)      = [];
    envelopeIntraCorrAR(removeDataIdx,:) = [];
    envelopeIntraCorrBR(removeDataIdx,:) = [];
    infraIntraCorrAR(removeDataIdx,:)    = [];
    infraIntraCorrBR(removeDataIdx,:)    = [];

    medPairCorrSuperR(removeDataIdx,:)      = [];
    medPairCorrMidR(removeDataIdx,:)        = [];
    medPairCorrDeepR(removeDataIdx,:)       = [];
    envelopePairCorrSuperR(removeDataIdx,:) = [];
    envelopePairCorrMidR(removeDataIdx,:)   = [];
    envelopePairCorrDeepR(removeDataIdx,:)  = [];
    infraPairCorrSuperR(removeDataIdx,:)    = [];
    infraPairCorrMidR(removeDataIdx,:)      = [];
    infraPairCorrDeepR(removeDataIdx,:)     = [];

    medPairASuperBMidR(removeDataIdx,:)  = [];
    medPairASuperBDeepR(removeDataIdx,:) = [];
    medPairAMidBSuperR(removeDataIdx,:)  = [];
    medPairAMidBDeepR(removeDataIdx,:)   = [];
    medPairADeepBSuperR(removeDataIdx,:) = [];
    medPairADeepBMidR(removeDataIdx,:)   = [];

    envelopeASuperBMidR(removeDataIdx,:)  = [];
    envelopeASuperBDeepR(removeDataIdx,:) = [];
    envelopeAMidBSuperR(removeDataIdx,:)  = [];
    envelopeAMidBDeepR(removeDataIdx,:)   = [];
    envelopeADeepBSuperR(removeDataIdx,:) = [];
    envelopeADeepBMidR(removeDataIdx,:)   = [];

    infraASuperBMidR(removeDataIdx,:)  = [];
    infraASuperBDeepR(removeDataIdx,:) = [];
    infraAMidBSuperR(removeDataIdx,:)  = [];
    infraAMidBDeepR(removeDataIdx,:)   = [];
    infraADeepBSuperR(removeDataIdx,:) = [];
    infraADeepBMidR(removeDataIdx,:)   = [];

    intraCorrAR(removeDataIdx,:) = [];    
    intraCorrBR(removeDataIdx,:) = [];  

    envelopeIntraAAllR(removeDataIdx,:) = []; 
    envelopeIntraBAllR(removeDataIdx,:) = []; 

    infraIntraAAllR(removeDataIdx,:) = [];
    infraIntraBAllR(removeDataIdx,:) = [];

    medPairCorrFreqTimeR(:,removeDataIdx)                = [];
    medPairCorrFreqTimeIntraAR(:,removeDataIdx)          = [];
    medPairCorrFreqTimeIntraBR(:,removeDataIdx)          = [];

    medPairCorrFreqPowerR(:,removeDataIdx)                = [];
    medPairCorrFreqPowerIntraAR(:,removeDataIdx)          = [];
    medPairCorrFreqPowerIntraBR(:,removeDataIdx)          = [];

    medPairCorrFreqInfraR(:,removeDataIdx)       = [];
    medPairCorrFreqInfraIntraAR(:,removeDataIdx) = [];
    medPairCorrFreqInfraIntraBR(:,removeDataIdx) = [];
    
freqCombNames = {'Theta-Alpha','Theta-Beta','Theta-Gamma','Theta-Spiking',...
    'Alpha-Beta','Alpha-Gamma','Alpha-Spiking','Beta-Gamma','Beta-Spiking','Gamma-Spiking'}; 


    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat'],'connValsR','distValsR',...
        'medPairCorrR','medCorrEnvelopeR','medCorrInfraSlowR','medIntraCorrAR','medIntraCorrBR','envelopeIntraCorrAR',...
        'envelopeIntraCorrBR','infraIntraCorrAR','infraIntraCorrBR','medPairCorrSuperR','medPairCorrMidR',...
        'medPairCorrDeepR','envelopePairCorrSuperR','envelopePairCorrMidR','envelopePairCorrDeepR',...
        'infraPairCorrSuperR','infraPairCorrMidR','infraPairCorrDeepR','medPairASuperBMidR','medPairASuperBDeepR',...
        'medPairAMidBSuperR','medPairAMidBDeepR','medPairADeepBSuperR','medPairADeepBMidR','envelopeASuperBMidR',...
        'envelopeASuperBDeepR','envelopeAMidBSuperR','envelopeAMidBDeepR','envelopeADeepBSuperR','envelopeADeepBMidR',...
        'infraASuperBMidR','infraASuperBDeepR','infraAMidBSuperR','infraAMidBDeepR','infraADeepBSuperR','infraADeepBMidR',...
        'removeDataIdx','intraCorrAR','intraCorrBR','envelopeIntraAAllR','envelopeIntraBAllR','infraIntraAAllR','infraIntraBAllR',...
        'medPairCorrFreqTimeR','medPairCorrFreqTimeIntraAR','medPairCorrFreqTimeIntraBR','medPairCorrFreqPowerR',...
        'medPairCorrFreqPowerIntraAR','medPairCorrFreqPowerIntraBR','medPairCorrFreqInfraR','medPairCorrFreqInfraIntraAR',...
        'medPairCorrFreqInfraIntraBR','freqCombNames','-append');
toc;

else
    disp('Loading saved variables...')
    clear allVars
    allVars = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\DualProbeVars.mat']);
    fieldNames = fieldnames(allVars);

    % Check if bad runs are completely removed from the processed data
    goodRunsR = reshape(goodRuns',[numel(goodRuns) 1]);

    % Remove recordings with a single channel
    singleChRow = cellfun(@(x) isscalar(x),allVars.intraCorrBR(:,1));

    goodRunsR(isnan(goodRunsR)) = 0;
    goodRunsR(allVars.removeDataIdx,:)  = [];

    for iL = 1:length(fieldNames)
        if ~(contains(fieldNames{iL},'Super') || contains(fieldNames{iL},'Mid') || contains(fieldNames{iL},'Deep'))
            allVars.(fieldNames{iL})(~goodRunsR,:) = [];
        else
            allVars.(fieldNames{iL})(~goodRunsR|singleChRow,:) = [];
        end
    end

end