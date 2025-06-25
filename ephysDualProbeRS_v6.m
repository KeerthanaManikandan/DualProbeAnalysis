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
monkeyName     = 'Whiskey';
hemisphere     = 'Left';
saveFigureFlag = 0;
fs             = 1e3; % Sampling frequency

gammaBand      = [30 90]; [zG,pG,kG] = butter(3,gammaBand./(fs/2),'bandpass'); [sosG,gG] = zp2sos(zG,pG,kG);
alphaBand      = [8 12];  [zA,pA,kA] = butter(3,alphaBand./(fs/2),'bandpass'); [sosA,gA] = zp2sos(zA,pA,kA);
betaBand       = [13 30]; [zB,pB,kB] = butter(3,betaBand./(fs/2),'bandpass');  [sosB,gB] = zp2sos(zB,pB,kB);
% probeBList     = [1:13 15:21 23:32];
chOutCortex    = 1:3;
chDeep         = 30:32;

disp(['Loading all required parameter names and values for ' monkeyName]);
[allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate,patchInfo,pairClass] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
clc; disp(['Loading all required parameter names and values for ' monkeyName ' ... Done']);
