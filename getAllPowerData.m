function [powA,powB,powEEG,specA,specB,specEEG,instPowA,instPowB,freqVals,freqValsSpec,timeValsSpec] = getAllPowerData(probeA,probeB,eeg,fs,monkeyName,expDate,fileNum)
% This function computes power spectrum density, spectrogram, instantaneous
% power (using hilbert transform) in the order of WB, alpha, beta, gamma 
% and stores/retrieves it. 

saveFolder = ['X:\Data\' monkeyName '_SqM\Left Hemisphere\' expDate '\Electrophysiology\All Power Data'];
if ~exist('saveFolder','dir'); [~,~] = mkdir(saveFolder); end

% Store/retrieve power for the two probes using pwelch
if ~exist([saveFolder '\allPowVals_' num2str(fileNum) '.mat'],'file')
    dispstat('','init');
    dispstat('Getting the power for the two probes... ');
    segLen = 500; winSize = 250;
    powA = zeros(120,size(probeA,2)); powB = zeros(120,size(probeB,2));

    for iCh = 1:size(probeA,2)
        [powA(:,iCh),freqVals] = pwelch(probeA(:,iCh),segLen, winSize,1:120,fs);
        if ~(iCh>size(probeB,2))
            [powB(:,iCh),freqVals] = pwelch(probeB(:,iCh),segLen,winSize,1:120,fs);
        end
        dispstat(['Getting the power for the two probes... ' num2str((iCh/size(probeA,2))*100) '% done']);
    end
    dispstat('Getting the power for the two probes... 100% done');

    % Get EEG power
    if ~isempty(eeg)
        powEEG = pwelch(eeg,segLen,winSize,1:120,fs);
    else
        powEEG = [];
    end

    % Get the spectrogram
    dispstat('','init');
    dispstat('Getting the spectrogram... ');
    params.Fs       = fs;
    params.fpass    = [1 120];
    params.pad      = -1;
    params.tapers   = [3 5];
    params.trialave = 0;

    if ~isempty(eeg);[specEEG,~, ~] = mtspecgramc(eeg,[5 2],params); else; specEEG = []; end

    [specA,~, ~]                      = mtspecgramc(probeA,[5 2],params);
    [specB,timeValsSpec,freqValsSpec] = mtspecgramc(probeB,[5 2],params);
    dispstat('Getting the spectrogram... Completed');
    %{
        % Plotting the spectrogram 
        figure;imagesc(time',freq',10.*log10(specgram'));
        set(gca,'YDir','normal')
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        colormap jet; shading interp; colorbar;
        caxis([ -20 70]);
    %}

    % Get the instantaneous power using Hilbert Transform for different
    % frequencies - can be updated to compute phase
    dispstat('','init');
    gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
    alphaBand = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
    betaBand  = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');  % Beta band filtering parameters
    bandLabels = {'Wideband';'Alpha band'; 'Beta band'; 'Gamma band'};

    instPowA  = cell(0); instPowB = cell(0);
    for iBand = 1:4
        switch iBand
            case 1 % WB
                xA = probeA; yA = probeB;

            case 2 % Alpha band
                xA = filtfilt(bA,aA,probeA);
                yA = filtfilt(bA,aA,probeB);

            case 3 % Beta band
                xA = filtfilt(bB,aB,probeA);
                yA = filtfilt(bB,aB,probeB);

            case 4 % Gamma band
                xA = filtfilt(bG,aG,probeA);
                yA = filtfilt(bG,aG,probeB);
        end
        dispstat(['Getting the instantaneous power for ' bandLabels{iBand} '... ']);
        instPowA{iBand,1}   = abs(hilbert(xA)).^2;
        instPowB{iBand,1}   = abs(hilbert(yA)).^2;
    end
    dispstat('Getting the instantaneous power for all bands ... Completed');
    
    % Save the power 
    disp('Saving the powers... ');
    save([saveFolder '\allPowVals_' num2str(fileNum) '.mat'],'powA','powB','powEEG','specA','specB','specEEG','instPowA',...
        'instPowB','freqVals','freqValsSpec','timeValsSpec');

else
    disp('Retrieving the powers... ');
    allVals      = load([saveFolder '\allPowVals_' num2str(fileNum) '.mat']); 
    powA         = allVals.powA; 
    powB         = allVals.powB;
    powEEG       = allVals.powEEG; 
    specA        = allVals.specA;
    specB        = allVals.specB;
    specEEG      = allVals.specEEG; 
    instPowA     = allVals.instPowA; 
    instPowB     = allVals.instPowB; 
    freqVals     = allVals.freqVals; 
    freqValsSpec = allVals.freqValsSpec;
    timeValsSpec = allVals.timeValsSpec; 

end
end