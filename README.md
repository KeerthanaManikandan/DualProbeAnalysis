# Read Me
This repository contains all codes necessary to process and store RS data collected from two linear electrode arrays simultaneously. Data are publicly available on the DANDI Archive (RRID:SCR_017571) at the following URL:

## Main scripts
1. `combinedAnalysisEphysDualProbeRS_vFinal.m`: Combines the pairwise correlations from two animals. 
2. 

## Main functions
1. `getDualProbeCorrelations.m`: Calculates the pairwise correltions between the electrodes for different timescales, frequencies, and compartments. Code also calculates cross-frequency correlations. This function calculates the correlations for a single animal 
2. `getPhaseAmpCoupling_v2.m`: Calculate the modulation index between the amplitude of gamma band and the phase of low frequencies (theta to beta)

## Dependent functions
1. `getAllDualProbeData.m`: Stores and retrieves LFP data, identifies bad channels and time segments. It also identifies the transition channels for both electrodes. 
2. `saveLFPDualProbe.m`: Saves the LFP data from the ns2 file. This also stores the EEG recorded from ripple. The neural signals are bandpass filtered between 6-250 Hz to remove artifacts due to respiration and heart rate. 
3. `getBadTimesAndChannels.m`: This function finds the bad time segments and channels; spectrograms are also calculated and saved to verify the removal. 

### Necessary folders
- 01_Sunin: To process imaging data.
- neuroshare: To process neural data recorded from Ripple
- nonlinear: To perform image registration.
- chronux_2_12: To calculate power spectral density, spectrograms. 

### Abbreviations
- rs-ISOI: Resting state intrinsic signal optical imaging
- LFP: Local Field Potential
- FC: Functional Connectivity 
- RS: Resting state
- ROI: Region of interest
- FOV: Field of view