function infraEphys = getInfraSlowPowerLFP(probeDat,bCoeff,aCoeff,chInCortex)
% Computes the infraslow fluctuations in power for the LFP - for different
% bands depending on the coefficients given as inputs to the signal 
% Updated on April 24, 2024 - KM 

% Filter the LFP recorded from channels in cortex and rectify it 
if length(chInCortex) == 2; chInCortex = chInCortex(1):chInCortex(2); end 

if ~isempty(bCoeff) && ~isempty(aCoeff)
    probeDat =  abs(single(filtfilt(bCoeff,aCoeff,double(probeDat(:,chInCortex)))));
else
    probeDat =  abs(probeDat(:,chInCortex));
end

% Get the envelope of the signal
envelopeDat = envelope(probeDat,5); %envelope(probeDat,5,'peak');

% Get filter parameters... 
clear z p k sos g enSize
fs      = 1e3;

% Bandpass - 0.01 Hz - 0.1 Hz
[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);
enSize  = size(envelopeDat); 

% Pad the envelope and filter the signals
envelopeFiltered = filtfilt(sos,g,double([envelopeDat; envelopeDat ;envelopeDat ]));
envelopeFiltered = envelopeFiltered(enSize(1)+1:(end-enSize(1)),:);
infraEphys       = single(downsample(envelopeFiltered,100));

% % Get the analytical version
% envelopeA = envelope(probeDat,5); 
% envelopeFilteredA = filtfilt(sos,g,double([envelopeA; envelopeA ;envelopeA ]));
% envelopeFilteredA = envelopeFilteredA(enSize(1)+1:(end-enSize(1)),:);
% infraEphysA       = single(downsample(envelopeFilteredA,100));

end