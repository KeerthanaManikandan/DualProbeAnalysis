function [modIdx,pj] = getPhaseAmpCoupling(lowFreq,highFreq,plotVar)
% Calculate the phase amplitude coupling.
% Based on Tort et al (2010), Journal of Neurophysiology 
% lowFreq: low frequency timecourse
% highFreq: high frequency timecourse
% Keerthana Manikandan
% October 13,2025

% Check if high and low frequency have the same number of channels
minChannels = min([size(highFreq,2) size(lowFreq,2)]); 

amplitude = abs(hilbert(highFreq(:,1:minChannels)));% Get amplitude of high frequency
phase     = angle(hilbert(lowFreq(:,1:minChannels))); % Get phase of low frequency

% Divide phase into bins
[~,edges,bin] = histcounts(rad2deg(phase),'BinMethod','integers','BinEdges',-180:20:180);

% Get the mean powers for each bin
n = size(edges,2)-1;

% Reshape 'bin' and 'amplitude' into column vectors
binCol = bin(:);
amplitudeCol = amplitude(:);

% Create a grouping index for each column
numColumns = size(bin, 2);
numRows = size(bin, 1);
groupIdx = repelem(1:numColumns, numRows)'; % repeats each column index for all its rows

% pj = mean(amplitude) at phase j
pj = accumarray([binCol, groupIdx], amplitudeCol, [], @mean);
pj = pj./sum(pj,1,'omitnan'); % normalized pac = pj/sum(pj for j= 1:n) % The sum should be equal to 1

% Calculate modulation index
hp = -sum(pj.*log(pj));
modIdx = (log(n)-hp)./log(n); 

% Plotting if desired
if ~exist('plotVar','var')|| isempty(plotVar)
    plotVar = 0; 
end

if plotVar
    % Plot the phase and amplitude in a polar plot
    figure; subplot(121); polarplot(median(phase,2,'omitnan'),median(amplitude,2,'omitnan'),'.');

    % Plot the distribution of means as a function of phase angle
    subplot(122); b =bar(edges(1:end-1)+180,median(pj,2,'omitnan'),'histc');ylim([0 0.1]);
    b.FaceColor = [0 0.4470 0.7410];
    xlim([-20 380]); xlabel('Phase'); ylabel('Mean PAC (normalized)');box off;
end

end