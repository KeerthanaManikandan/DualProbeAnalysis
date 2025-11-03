function [modIdx,pj] = getPhaseAmpCoupling(phaseVal,amplitudeVal,plotVar)
% Calculate the phase amplitude coupling.
% Based on Tort et al (2010), Journal of Neurophysiology 
% phaseVal: Low frequency phase
% amplitudeVal: High frequency amplitude
% Keerthana Manikandan
% October 13,2025

if nargin<3
    plotVar = 0; 
end 

% Divide phase into bins
[~,edges,bin] = histcounts(rad2deg(phaseVal),'BinEdges',-180:20:180);

% Get the mean powers for each bin
n = size(edges,2)-1;

% Reshape 'bin' and 'amplitude' into column vectors
binCol       = bin(:);
amplitudeCol = amplitudeVal(:);

% Create a grouping index for each column
numColumns = size(bin, 2);
numRows = size(bin, 1);
groupIdx = repelem(1:numColumns, numRows)'; % repeats each column index for all its rows

% pj = mean(amplitude) at phase j
pj = accumarray([binCol, groupIdx], amplitudeCol, [], @mean);
pj = pj./sum(pj,1,'omitnan'); % normalized pac = pj/sum(pj for j= 1:n) % The sum should be equal to 1

% Calculate modulation index
modIdx = (log(n)-(-sum(pj.*log(pj))))./log(n); 

if plotVar
    % Plot the phase and amplitude in a polar plot
    figure; subplot(121); polarplot(median(phaseVal,2,'omitnan'),median(amplitudeVal,2,'omitnan'),'.');

    % Plot the distribution of means as a function of phase angle
    subplot(122); b =bar(edges(1:end-1),median(pj,2,'omitnan'),'histc');ylim([0 0.1]);
    b.FaceColor = [0 0.4470 0.7410];
    xlim([-180 180]); xlabel('Phase'); ylabel('Mean PAC (normalized)');box off;
end

end