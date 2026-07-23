% Convert ephys data to NWB format for DANDI
% This code is for the Dual Probe paper
% July 22, 2026 - KM (with Claude's help)

% Refer here: https://docs.dandiarchive.org/getting-started/data-standards/#neurodata-without-borders-nwb
% Install MATNWB by cloning from GitHub : git clone https://github.com/NeurodataWithoutBorders/matnwb.git 
% Then add matnwb to path and run this function generateCore() to
% initialize NWB
% Refer: https://matnwb.readthedocs.io/en/latest/pages/tutorials/ecephys.html

% Set paths
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\DualProbe\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\DualProbe\neuroshare']));
addpath(genpath([commonDir '\Codes\DualProbe']));
addpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\DualProbe\chronux_2_12\spectral_analysis\continuous\dupes']));


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
            NaN NaN NaN NaN 0 0 0 0 0 0 NaN NaN NaN NaN NaN NaN; ...
            NaN NaN NaN NaN 0 0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN; ...
            0 1 1 0 1 0 1 1 1 NaN NaN NaN NaN NaN NaN NaN; ...
            1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
        ageVal = 'P15Y';
        monkeyNameShort = 'CS';

    case 2
        monkeyName = 'Whiskey';
        goodRuns = [1 1 1 1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; ...
            1 1 1 0 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN; ...
            1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
            1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 NaN; ...
            1 1 1 1 1 1 1 1 1 1 1 NaN NaN NaN NaN NaN;...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 NaN NaN]';
         ageVal = 'P6Y';
         monkeyNameShort = 'W';
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

%% Saving data in NWB format to share in DANDI database
clc;

for iDate = 1:size(allDates,1)
    clear expDate datFileNum saveFolder
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    dataDir = ['C:\Users\kem294\Documents\Data\DANDI\001881\sub-' monkeyNameShort];
    if ~exist(dataDir,'dir'); [~,~] = mkdir(dataDir); end


    for iRun = 1:length(datFileNum)
        clear fileNum recording
        fileNum = datFileNum(iRun);


        if ~goodRuns(fileNum,iDate) || isnan(goodRuns(fileNum,iDate)); continue; end


        % delete existing files if present
        if exist(fullfile(dataDir,[monkeyNameShort '_' expDate '_' fileNum '_nwb.mat']), 'file')
            delete(fullfile(dataDir,[monkeyNameShort '_' expDate '_' fileNum '_nwb.mat']));
        end
        if exist(fullfile(dataDir,[monkeyNameShort '_' expDate '_' fileNum '.nwb']), 'file')
            delete(fullfile(dataDir,[monkeyNameShort '_' expDate '_' fileNum '.nwb']));
        end

        disp(['Saving data for ' monkeyNameShort ' Exp: ' expDate ' run: ' num2str(fileNum)]);

        % Build recording struct
        % ephys data — stored as [time x channels]
        recording.lfp1           = single(allProbeData{fileNum,iDate}.probe1Ch)';
        recording.lfp2           = single(allProbeData{fileNum,iDate}.probe2Ch)';
        recording.raw1           = single(allProbeData{fileNum,iDate}.raw1Ch)';
        recording.raw2           = single(allProbeData{fileNum,iDate}.raw2Ch)';
        recording.eeg            = single(allProbeData{fileNum,iDate}.eegCh)';
        recording.fs             = 1e3;

        % channel info
        recording.channels1      = single(1:size(allProbeData{fileNum,iDate}.probe1Ch,2));
        recording.channels2      = single(1:size(allProbeData{fileNum,iDate}.probe2Ch,2));

        recording.bad_channels1  = single(badElecA{fileNum,iDate});
        recording.bad_channels2  = single(badElecA{fileNum,iDate});

        recording.bad_times      = single(allBadTimes{fileNum,iDate});
        if isempty(recording.bad_times); recording.bad_times = 0; end

        recording.cortex_channels1 = single(estChInCortexA{iDate}(iRun,1):estChInCortexA{iDate}(iRun,2));
        recording.cortex_channels2 = single(estChInCortexB{iDate}(iRun,1):estChInCortexB{iDate}(iRun,2));

        % inter-electrode measurements
        recording.distance         = distSites{iDate}(fileNum);      % scalar, distance in mm between probes
        recording.connectivity     = connValsAll(fileNum,iDate);

        % session info
        expDateNew = datetime(expDate, 'InputFormat', 'MM_dd_yyyy');
        expDateNew = char(datetime(expDateNew, 'Format', 'yyyyMMdd'));

        recording.subject        = monkeyNameShort;
        recording.session_id     = [expDateNew fileNum];
        recording.date           = expDateNew;
        recording.brain_region1  = pairClass{iDate}(iRun,1);
        recording.brain_region2  = pairClass{iDate}(iRun,2);

        % Create NwbFile
        nwbfile = NwbFile(...
            'identifier',               recording.session_id,...
            'session_description',      'Simultaneous RS electrophysiology using two linear electrodes',...
            'session_start_time',       datetime([recording.date ' 00:00:00'], 'InputFormat', 'yyyyMMdd HH:mm:ss', 'TimeZone', 'local'),...
            'general_institution',      'University of Pittsburgh',...
            'general_lab',              'Gharbawie Lab',...
            'general_protocol',         'R01NS143962',...
            'general_source_script',    'https://github.com/KeerthanaManikandan/DualProbeAnalysis',...
            'general_source_script_file_name', 'nwbDataConversion_DualProbe.m');

        % Subject metadata
        nwbfile.general_subject = types.core.Subject(...
            'subject_id',   recording.subject,...
            'species',      'http://purl.obolibrary.org/obo/NCBITaxon_27679',...
            'sex',          'M',...
            'age',          ageVal);

        % Devices — one per probe
        device1 = types.core.Device(...
            'description',  'A1x32-15 mm',...
            'manufacturer', 'Neuronexus');
        nwbfile.general_devices.set('probe1', device1);

        device2 = types.core.Device(...
            'description',  'A1x32-15 mm',...
            'manufacturer', 'Neuronexus');
        nwbfile.general_devices.set('probe2', device2);

        % Electrode groups — one per probe
        egroup1 = types.core.ElectrodeGroup(...
            'description',  'Electrode 1',...
            'location',     recording.brain_region1,...
            'device',       types.untyped.SoftLink(device1));
        nwbfile.general_extracellular_ephys.set('shank1', egroup1);

        egroup2 = types.core.ElectrodeGroup(...
            'description',  'Electrode 2',...
            'location',     recording.brain_region2,...
            'device',       types.untyped.SoftLink(device2));
        nwbfile.general_extracellular_ephys.set('shank2', egroup2);

        % Electrodes table — one row per channel across both probes
        nchans1 = length(recording.channels1);
        nchans2 = length(recording.channels2);

        location_col1   = repmat({recording.brain_region1}, nchans1, 1);
        bad_col1        = ismember(recording.channels1, recording.bad_channels1)';
        group_col1      = repmat(types.untyped.ObjectView(egroup1), nchans1, 1);
        group_name_col1 = repmat({'shank1'}, nchans1, 1);

        location_col2   = repmat({recording.brain_region2}, nchans2, 1);
        bad_col2        = ismember(recording.channels2, recording.bad_channels2)';
        group_col2      = repmat(types.untyped.ObjectView(egroup2), nchans2, 1);
        group_name_col2 = repmat({'shank2'}, nchans2, 1);

        nwbfile.general_extracellular_ephys_electrodes = types.core.ElectrodesTable(...
            'description',  'electrode metadata',...
            'colnames',     {'location', 'bad', 'group', 'group_name'},...
            'id',           types.hdmf_common.ElementIdentifiers('data', int64(0:nchans1+nchans2-1)'),...
            'location',     types.hdmf_common.VectorData('data', [location_col1;   location_col2],   'description', 'brain region'),...
            'bad',          types.hdmf_common.VectorData('data', [bad_col1;        bad_col2],        'description', 'bad channel flag'),...
            'group',        types.hdmf_common.VectorData('data', [group_col1;      group_col2],      'description', 'electrode group'),...
            'group_name',   types.hdmf_common.VectorData('data', [group_name_col1; group_name_col2], 'description', 'electrode group name'));

        % Electrode table regions
        % probe 1: rows 0 to nchans1-1
        electrode_table_region1 = types.hdmf_common.DynamicTableRegion(...
            'table',        types.untyped.ObjectView(nwbfile.general_extracellular_ephys_electrodes),...
            'data',         int32(0:nchans1-1)',...
            'description',  'probe 1 electrodes');

        % probe 2: rows nchans1 to nchans1+nchans2-1
        electrode_table_region2 = types.hdmf_common.DynamicTableRegion(...
            'table',        types.untyped.ObjectView(nwbfile.general_extracellular_ephys_electrodes),...
            'data',         int32(nchans1:nchans1+nchans2-1)',...
            'description',  'probe 2 electrodes');

        % EEG uses row 0 (first electrode as reference)
        eeg_table_region = types.hdmf_common.DynamicTableRegion(...
            'table',        types.untyped.ObjectView(nwbfile.general_extracellular_ephys_electrodes),...
            'data',         int32(0),...
            'description',  'EEG electrode');

        % Convert data from microvolts to volts
        recording.raw1 = double(recording.raw1) / 1e6;
        recording.raw2 = double(recording.raw2) / 1e6;
        recording.lfp1 = double(recording.lfp1) / 1e6;
        recording.lfp2 = double(recording.lfp2) / 1e6;
        recording.eeg  = double(recording.eeg)  / 1e6;

        % Store raw and LFP for probe 1 and probe 2 in processing/ecephys
        raw_es1 = types.core.ElectricalSeries(...
            'data',               recording.raw1,...
            'starting_time',      0.0,...
            'starting_time_rate', double(recording.fs),...
            'electrodes',         electrode_table_region1,...
            'description',        'probe 1 broadband signal downsampled to 1kHz');

        raw_es2 = types.core.ElectricalSeries(...
            'data',               recording.raw2,...
            'starting_time',      0.0,...
            'starting_time_rate', double(recording.fs),...
            'electrodes',         electrode_table_region2,...
            'description',        'probe 2 broadband signal downsampled to 1kHz');

        lfp_es1 = types.core.ElectricalSeries(...
            'data',               recording.lfp1,...
            'starting_time',      0.0,...
            'starting_time_rate', double(recording.fs),...
            'electrodes',         electrode_table_region1,...
            'description',        'probe 1 LFP 6-250 Hz bandstop at 60 Hz');

        lfp_es2 = types.core.ElectricalSeries(...
            'data',               recording.lfp2,...
            'starting_time',      0.0,...
            'starting_time_rate', double(recording.fs),...
            'electrodes',         electrode_table_region2,...
            'description',        'probe 2 LFP 6-250 Hz bandstop at 60 Hz');

        lfp_obj = types.core.LFP(...
            'probe1', lfp_es1,...
            'probe2', lfp_es2);

        ecephys = types.core.ProcessingModule(...
            'description',  'ecephys',...
            'LFP',          lfp_obj,...
            'raw_probe1',   raw_es1,...
            'raw_probe2',   raw_es2);

        nwbfile.processing.set('ecephys', ecephys);

        % Store EEG
        if ~isempty(recording.eeg) && any(recording.eeg(:) ~= 0)
            eeg_es = types.core.ElectricalSeries(...
                'data',               recording.eeg,...
                'starting_time',      0.0,...
                'starting_time_rate', double(recording.fs),...
                'electrodes',         eeg_table_region,...
                'description',        'Intracranial EEG');
            nwbfile.acquisition.set('EEG', eeg_es);
        end

        % eeg_es = types.core.ElectricalSeries(...
        %     'data',               recording.eeg,...
        %     'starting_time',      0.0,...
        %     'starting_time_rate', double(recording.fs),...
        %     'electrodes',         eeg_table_region,...
        %     'description',        'Intracranial EEG');
        % nwbfile.acquisition.set('EEG', eeg_es);

        % Store bad times
        bad_ts = types.core.TimeSeries(...
            'data',         recording.bad_times,...
            'timestamps',   recording.bad_times / recording.fs,...
            'description',  'indices of bad time samples',...
            'data_unit',    'samples');
        nwbfile.acquisition.set('bad_times', bad_ts);

        % Store inter-electrode distance and connectivity as TimeSeries
        distance_ts = types.core.TimeSeries(...
            'data',         recording.distance,...
            'timestamps',   0.0,...
            'description',  'distance between probe 1 and probe 2 in mm',...
            'data_unit',    'mm');
        nwbfile.acquisition.set('inter_electrode_distance', distance_ts);

        connectivity_ts = types.core.TimeSeries(...
            'data',         recording.connectivity,...
            'timestamps',   0.0,...
            'description',  'connectivity values between probe 1 and probe 2',...
            'data_unit',    'a.u.');
        nwbfile.acquisition.set('functional_connectivity', connectivity_ts);

        % Export NWB file
        nwbExport(nwbfile, fullfile(dataDir, ['sub-' monkeyNameShort '_ses-' expDateNew 'run' num2str(fileNum) '_ecephys.nwb']));
    end
end

%% How to read the NWB data
nwb = nwbRead(fullfile(dataDir, ['sub-' monkeyNameShort '_ses-' expDateNew 'run' num2str(fileNum) '_ecephys.nwb']));

% LFP
lfp_probe1 = nwb.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('probe1').data;
lfp_probe2 = nwb.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('probe2').data;

% raw
raw_probe1 = nwb.processing.get('ecephys').nwbdatainterface.get('raw_probe1').data;
raw_probe2 = nwb.processing.get('ecephys').nwbdatainterface.get('raw_probe2').data;

% EEG
eeg = nwb.acquisition.get('EEG').data;

% bad times
bad_times = nwb.acquisition.get('bad_times').data;

% distance and connectivity
distance     = nwb.acquisition.get('inter_electrode_distance').data;
fc = nwb.acquisition.get('functional_connectivity').data;

% electrodes table
elecTable = nwb.general_extracellular_ephys_electrodes;

