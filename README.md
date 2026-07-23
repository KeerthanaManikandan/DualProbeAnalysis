# Read Me
This repository contains all codes necessary to process and store resting state electrophysiology collected from two linear electrode arrays simultaneously. Data are publicly available on the DANDI Archive (RRID:SCR_017571) at the following URL: https://dandiarchive.org/dandiset/001881

## Main scripts
1. `combinedAnalysisEphysDualProbeRS_vFinal.m`: Combines the pairwise correlations from two animals and generates the necessary figures in the paper
2. `ephysDualProbeRS_vFinal.m`: Stores or retrieves the LFP and spiking data for each recording for a single animal. This script also calculates the pairwise correlations and coherence between neural data from the two linear electrodes. 

## Main functions
1. `getAllDualProbeData.m`: Stores and retrieves LFP data, identifies bad channels and time segments. It also identifies the transition channels for both electrodes. 
2. `saveLFPDualProbe.m`: Saves the LFP data from the ns2 file. This also stores the EEG recorded from ripple. The neural signals are bandpass filtered between 6-250 Hz to remove artifacts due to respiration and heart rate. 
3. `getBadTimesAndChannels.m`: This function finds the bad time segments and channels; spectrograms are also calculated and saved to verify the removal. 
4. `getDualProbeCorrelations.m`: Calculates the pairwise correltions between the electrodes for different timescales, frequencies, and compartments. Code also calculates cross-frequency correlations. This function calculates the correlations for a single animal. 
5. `getMonkeyParamsDualProbeEphys.m`: Obtains the experiment dates, recording numbers, electrode labels and so on for a single animal. 
6. `getDistanceAndConnValsDualProbe.m`: This function outputs the distance between the electrodes (in mm) and also calculates the functional connectivity value with respect to the reference electrode. 
7. `getRSConnectivityMaps.m`: Calculates the FC map for a given seed. 
8. `nwbDataConversion.m`: Converts the LFP into format that is compatible with DANDI database.


### Necessary folders
- 01_Sunin: To process imaging data.
- neuroshare: To process neural data recorded from Ripple
- nonlinear: To perform image registration.
- chronux_2_12: To calculate power spectral density, spectrograms and coherence
- matnwb: To convert electrophysiological data into Neurodata without borders format

### Abbreviations
- rs-ISOI: Resting state intrinsic signal optical imaging
- LFP: Local Field Potential
- FC: Functional Connectivity 
- RS: Resting state
- ROI: Region of interest
- FOV: Field of view