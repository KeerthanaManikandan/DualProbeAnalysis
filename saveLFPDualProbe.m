function [probe1Ch,probe2Ch,raw1Ch,raw2Ch,eegCh,scalpEEGCh] = saveLFPDualProbe(monkeyName,expDate,fileNum,datFileName,serverPath,fs)

% This function saves the LFP data (1-250 Hz) for the two probes for a
% specific file in a particular experiment. The EEG and/or scalp EEG is
% also stored. The LFP is bandpassed from 6-250 Hz to remove artifacts
% caused due to respiration and heart rate

clc;clear probe1Ch probe2Ch eegCh scalpEEGCh probeCh rawCh raw1Ch raw2Ch
probe1Ch = []; probe2Ch = []; eegCh = []; scalpEEGCh = [];
raw1Ch = []; raw2Ch = [];

disp(['Obtaining and processing data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
% datName = [serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum) ];

if exist([serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum) '.nev' ],'file')
    datName = [serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum)];

elseif exist([serverPath expDate '\Electrophysiology\run0' num2str(fileNum) '\' datFileName num2str(fileNum) '.nev'],'file')
    datName = [serverPath expDate '\Electrophysiology\run0' num2str(fileNum) '\' datFileName num2str(fileNum)];
end

% Access and open data file
try
    [nsResult,hFile] = ns_OpenFile(datName);
catch
    disp('Data did not load.... Moving on to the next recording');
    return;
end

if ~strcmp(nsResult,'ns_OK')
    datName = [serverPath expDate '\Electrophysiology\run0' num2str(fileNum-1) '\' datFileName num2str(fileNum)];
    [nsResult,hFile] = ns_OpenFile(datName);

    if ~strcmp(nsResult,'ns_OK')
        datName = [serverPath expDate '\Electrophysiology\run' num2str(fileNum-1) '\' datFileName num2str(fileNum)];
        [nsResult,hFile] = ns_OpenFile(datName);
    end

    if ~strcmp(nsResult,'ns_OK')
        disp('Data file did not open! - going to the next datafile');
        return;
    end
end

% Get file information
[nsResult2, fileInfo] = ns_GetFileInfo(hFile);
if ~strcmp(nsResult2,'ns_OK')
    disp('Data file information did not load!');
    return;
end

% Get entity information
for iEntity = 1: fileInfo.EntityCount
    [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity); %#ok<AGROW> 
end

% Sort the entities whether they contain events, neural data, segment data
% Get the label  and list of all the channels storing LFP data
allList   = find([entityInfo.EntityType] == 2);
allLabel  = {entityInfo(allList).EntityLabel}; 

% Get indices of the LFP only
if strcmp(allLabel{1}(1:3),'lfp')
    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(1:3),'lfp'),allLabel,'un',0)); % 32-channel electrode

elseif strcmp(allLabel{1}(end-2:end),'lfp')
    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(end-2:end),'lfp'),allLabel,'un',0)); % 32-channel electrode

elseif strcmp(allLabel{1}(1:6),'analog')
    lfpIdx = cell2mat(cellfun(@(x) strcmp(x,'analog 2'),allLabel,'un',0)); % Single electrode
end

lfpList  = allList(lfpIdx);
lfpLabel = allLabel(lfpIdx);

if strcmp(lfpLabel{1},lfpLabel{2}) % To remove 30kHz sampled data
    lfpList(2) = [];
    lfpLabel(2) = [];
end

if length(lfpList)~= 1
    rawList = allList(~lfpIdx);
    rawLabel = allLabel(~lfpIdx);
else
    rawList = allList(end); % Single channel electrode condition 
    rawLabel = allLabel(end);
end

% Remove EEG label from raw data list
if find(cell2mat((cellfun(@(x)(strcmp(x(1:5),'analo')),rawLabel,'un',0)))) & length(rawList)~= 1 % To include for single electrode
    loc = find(cell2mat((cellfun(@(x)(strcmp(x(1:5),'analo')),rawLabel,'un',0)))); 
    rawList(loc)  = []; 
    rawLabel(loc) = []; 
end 

clear uniqueCh

% Find if there are two sets of 32 channel electrodes or not
uniqueCh = unique(cellfun(@(x) x(1),lfpLabel,'UniformOutput',0)');
allSets  = cellfun(@(x) x(1),lfpLabel,'UniformOutput',0)';
namesNotListedFlag = 0;

% If the two sets of electrodes are not properly labeled
if sum(lfpIdx) == 0 && (~strcmp(expDate,'08_08_2022')) % Removing the condition where the dual probe was done with one 32 Ch and one single electrode
    lfpList   = find([entityInfo.EntityType] == 2);
    lfpLabel  = {entityInfo(lfpList).EntityLabel}; % Gets the label of all the channels in lfpList
    lfpIdx    = cell2mat(cellfun(@(x) strcmp(x(1:3),'lfp'),lfpLabel,'un',0));

    lfpList(~lfpIdx)  = [];
    lfpLabel(~lfpIdx) = [];

    uniqueCh = {'A';'B'};
    namesNotListedFlag = 1;

elseif isscalar(uniqueCh) % If only one label was identified
    uniqueCh = {'A';'B'};
    namesNotListedFlag = 1;   
end

% Get the LFP for the probes recorded
disp('Getting the LFP for the two probes and filtering them ... ');

% Set filtering parameters
clear b a bS aS
[b,a] = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
[bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz
[bH,aH] = butter(3,250./(30e3/2),'high'); % High pass filtering >250 Hz

for iSet = 1:size(uniqueCh,1)
    clear setIdx channels chNum lfpListChWise probe rawCh
    if ~namesNotListedFlag
        setIdx = find(strcmp(allSets,uniqueCh(iSet)));
    else
        if iSet == 1
            setIdx = 1:32;
        else
            setIdx = 33:64;
        end
    end

    channels      = lfpLabel(setIdx);
    lfpListChWise = lfpList(setIdx);
    rawListChWise = rawList(setIdx);

    if ~namesNotListedFlag
        chNum = str2num(cell2mat(cellfun(@(x) x(2:3),channels','un',0))); %#ok<ST2NM> 
    else
        chNum = str2num(cell2mat(cellfun(@(x) x(end-1:end),channels','un',0))); %#ok<ST2NM> 
    end

    [~,elecID] = sort(chNum);  % Rearranging the channels based on the channel map

    % Get the LFP for all the channels in sorted order
    for iElec = 1:length(channels)
        clear elecEntityID lfpEntityID lfpCount
        if ~isempty(elecID)
            lfpEntityID   = lfpListChWise(elecID(iElec));
            rawEntityID   = rawListChWise(elecID(iElec));
        else
            lfpEntityID   = lfpListChWise(iElec);
            rawEntityID   = rawListChWise(iElec);
        end

        lfpCount      = entityInfo(lfpEntityID).ItemCount;
        rawCount      = entityInfo(rawEntityID).ItemCount;

        % Get LFP data
        [~, ~, probeCh(:,iElec)] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
        [~, ~, rawCh(:,iElec)]   = ns_GetAnalogData(hFile,rawEntityID,1,rawCount);
    end

    probeCh = filtfilt(b,a,probeCh); % Bandpass filtering across 1-250 Hz
    probeCh = filtfilt(bS,aS,probeCh); % Bandstop filtering between 57-62 Hz
    rawCh   = single(downsample(filtfilt(bH,aH,rawCh),30)); % Highpass frequencies beyond 250 Hz

    switch iSet
        case 1
            probe1Ch = probeCh;
            raw1Ch   = rawCh;
        case 2
            probe2Ch = probeCh;
            raw2Ch   = rawCh; 
    end
end

% Getting the EEG and/or single probe
clear elecList elecID
elecList   = find([entityInfo.EntityType] == 2);
elecLabel  = {entityInfo(elecList).EntityLabel}; % Gets the label of all the channels in lfpList
elecIdx    = cell2mat(cellfun(@(x) strcmp(x(1:end-2),'analog'),elecLabel,'un',0));
elecList(~elecIdx) = [];
elecLabel(~elecIdx) = [];

if ~isempty(elecList)
    for iC = 1:length(elecList)
        clear elecCount analogInfo2 fs_ns5 
        elecCount      = entityInfo(elecList(iC)).ItemCount;
        [~, analogInfo2] = ns_GetAnalogInfo(hFile, elecList(iC));
        fs_ns5  = analogInfo2.SampleRate;
        if strcmp(elecLabel{iC},'analog 2') % analog 2 has the single probe
            clear probe2Temp
            [~, ~, probe2Temp] = ns_GetAnalogData(hFile,elecList(iC),1,elecCount);

            if ~strcmp(expDate,'10_17_2022')
                % Downsample the analog 2 channel from 30kHz to
                % 1kHz and bandpass filter across 1- 250 Hz
                probe2Ch = downsample(probe2Temp,fs_ns5/fs);
                probe2Ch= filtfilt(b,a,probe2Ch);
                probe2Ch = filtfilt(bS,aS,probe2Ch); % Bandstop filtering between 57-62 Hz

                % Get spiking frequencies
                raw2Ch = single(downsample(filtfilt(bH,aH,probe2Temp),30)); 

            else % Scalp EEG is collected in analog 2 for Whiskey on 10_17_2022
                scalpEEGCh = downsample(probe2Temp,fs_ns5/fs);
                scalpEEGCh = filtfilt(b,a,scalpEEGCh);
            end

        elseif strcmp(elecLabel{iC},'analog 1') % analog 1 has EEG
            [~, ~, eegCh] = ns_GetAnalogData(hFile,elecList(iC),1,elecCount);
            % Downsample if EEG is stored at a sampling
            % rate of 30kHz
            if fs_ns5 == 30e3; eegCh = downsample(eegCh,fs_ns5/fs); end
            eegCh = filtfilt(b,a,eegCh);
        end
    end
end

% Filter the LFP to remove 1-5 Hz.
if ~exist('fs','var'); fs = 1e3; end

% Filter the LFP to remove 1-5 Hz
[bL,aL] = butter(3,([6 250]./(fs/2)),'bandpass');
probe1Ch = single(filtfilt(bL,aL,double(probe1Ch)));
probe2Ch = single(filtfilt(bL,aL,double(probe2Ch)));

% Remove 60 Hz power line noise from EEG
eegCh = single(filtfilt(bS,aS,eegCh)); % Bandstop filtering between 57-62 Hz

if ~isempty(scalpEEGCh); scalpEEGCh = single(filtfilt(bS,aS,scalpEEGCh)); end 
end