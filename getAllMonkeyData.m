% This code retrieves data from all animals (here, Charlie and Whiskey) and
% further analysis is done on the data
% March 03, 2023 - Keerthana Manikandan
% Begin initialization
clear;clc;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir));
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));

hemisphere = 'Left'; 
allSqM     = {'CharlieSheen'; 'Whiskey'};
segLen     = 500;
winSize    = 250;

for iM = 1:size(allSqM,1)
    clear monkeyName allDates datFileNumAll serverPath refDate refImageName datFileNumAll chInCortexProbeA ...
        chInCortexProbeB chInCortexProbeAUpdated chInCortexProbeBUpdated
    monkeyName = allSqM{iM};
    % 1. Get monkey params
    [allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    chInCortexProbeAUpdated,chInCortexProbeBUpdated] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
    
   % 2. Store/Retrieve the LFP (1-250 Hz) and also removing 1-5 Hz
   disp(['Storing/Retrieving LFP data for ' monkeyName]);
   [probe1{iM},probe2{iM},eeg{iM},scalpEEG{iM}] = saveLFPDualProbe(monkeyName,hemisphere,allDates,datFileNameAll,datFileNumAll,serverPath);
   disp('Data Stored/Retrieved');

   % 3.Get/retrieve the distance and the connectivity values for the sites
   clear distSites connSites greenMapRef
   [distSites{iM},connVals{iM},refSites{iM},movSites{iM},greenMapRef{iM}] = getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName);

   
end 


tic;
dispstat('','init');
dispstat('Getting power for EEGs... ');
count = 1; 
for iDate = 1:size(eegWhiskey,2)
    for iFile = 1:size(eegWhiskey,1)
        if isempty(eegWhiskey{iFile,iDate})
            count = count+1; 
            continue; 
        else
            [powEEGWhiskey{iFile,iDate},f] = pwelch(eegWhiskey{iFile,iDate},segLen, winSize,1:120,fs);
            count = count+ 1;         
        end 
           dispstat(['Getting power for EEGs...' num2str((count/numel(eegWhiskey))*100) '% Complete']);
    end 
end 
dispstat('Getting power for EEGs... 100% Complete');
toc;
%%
figure; 
for iDate = 1:size(eegCharlie,2)
    for iFile = 1:size(eegCharlie,1)
        if isempty(powEEGCharlie{iFile,iDate})
            continue; 
        else
            plot(f,10.*log10(powEEGCharlie{iFile,iDate}));  hold on; 
        end 
    end
end 