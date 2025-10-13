function [modIdx,normP] = getPhaseAmpCoupling(lowFreq,highFreq,plotVar)
% Calculate the phase amplitude coupling based on Tort et al 2010

amplitude = abs(hilbert(highFreq)); % Get amplitude of high frequency
phase     = angle(hilbert(lowFreq)); % Get phase of low frequency

% Divide phase in 20 degree bins
[~,edges,bin] = histcounts(rad2deg(phase),'BinMethod','integers','BinEdges',-180:20:180);

% Get the mean powers for each bin
n = size(bin,2);
for iM = 1:n
    pj(:,iM) = accumarray(bin(:,iM),amplitude(:,iM),[],@mean); % pj = mean(amplitude) at phase j
end
normP = pj./sum(pj,2,'omitnan'); % normalized pac = pj/sum(pj for j= 1:n)

% Calculate modulation index
hp = -sum(pj.*log2(pj));
modIdx = (log2(n)-hp)./log2(n); 

% Plotting if desired
if ~exist('plotVar','var')|| isempty(plotVar)
    plotVar = 0; 
end

if plotVar
    % Plot the phase and amplitude in a polar plot
    figure; subplot(121); polarplot(median(phase,2,'omitnan'),median(amplitude,2,'omitnan'),'.');

    % Plot the distribution of means as a function of phase angle
    subplot(122); b =bar(edges(1:end-1)+180,median(normP,2,'omitnan'),'histc');ylim([0 0.1]);
    b.FaceColor = [0 0.4470 0.7410];
    xlim([-20 380]); xlabel('Phase'); ylabel('Mean PAC (normalized)');box off;
end

end