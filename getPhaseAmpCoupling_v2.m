function modIdx = getPhaseAmpCoupling_v2(phaseVal,amplitudeVal,winStart,winEnd,winSize,ampIdx,phaseIdx)

% Calculate the phase amplitude coupling.
% Based on Tort et al (2010), Journal of Neurophysiology 
% Keerthana Manikandan
% October 13,2025
% Updated on Jan 12, 2026 

% Initialize variables
nLow     = size(phaseVal,1);     % # Low frequency windows 
nHigh    = size(amplitudeVal,1); % # High frequency windows
nWin     = numel(winStart);    % Windows to iterate
nChan    = size(ampIdx,1);     % Channel combinations

nBins    = 18;
binEdges = -180:20:180;
logN     = single(log(nBins)); 
groupIdx = single(repelem(1:nChan,winSize)');
phaseVal   = rad2deg(phaseVal);

modIdx   = NaN(nWin,nHigh,nLow,nChan,'single'); % Initialization

% Progress tracking - Uncomment the code below if you want to track code
% progress within parfor loop
% D = parallel.pool.DataQueue;
% progress = 0;
% startTime = tic;
% nChars    = 1;
% afterEach(D, @updateProgress);
% 
%     function updateProgress(~)
%         progress = progress + 1;
%         elapsed  = toc(startTime);
%         pct = 100 * progress / nWin;
%         avgTime = elapsed / progress;
%         remaining = (nWin - progress) * avgTime;
% 
%         % Clear previous line and print update
%         fprintf(repmat('\b', 1, nChars));
%         nChars = fprintf('\rProgress: %d/%d (%.1f%%) | Elapsed: %.1fm | Remaining: %.1fm | Avg: %.1fs/win', ...
%             progress, nWin, pct, elapsed/60, remaining/60, avgTime);
%     end

parfor iWin = 1:nWin 
    amplitudeTemp = squeeze(amplitudeVal(:,winStart(iWin):winEnd(iWin),ampIdx));  % Obtain amplitude
    phaseTemp     = squeeze(phaseVal(:,winStart(iWin):winEnd(iWin),phaseIdx)); % Obtain phase %#ok<*PFBNS>

    [~,~,binTemp]     = histcounts(phaseTemp,'BinEdges',binEdges); binTemp = single(binTemp); % Obtain histogram 

    modIdxTemp    = NaN(nHigh, nLow, nChan, 'single'); 
   
    for iHigh = 1:nHigh
        % Get the amplitude for high frequency
        ampTempHigh      = squeeze(amplitudeTemp(iHigh,:,:));
        amplitudeCol = ampTempHigh(:);     

        for iLow = 1:nLow
            % Get the phase for low frequency 
            binTempLow = squeeze(binTemp(iLow,:,:));
            binTempCol = binTempLow(:);         
            
            pj = accumarray([binTempCol, groupIdx], amplitudeCol, [], @mean);
            pj = pj./sum(pj,1,'omitnan'); % The sum should be equal to 1

            % Calculate the KL divergence between amplitude of high
            % frequency and the phase of low frequency
            modIdxTemp(iHigh, iLow, :)  = (logN-(-sum(pj.*log(pj))))./logN; 
        end
    end

    % Assign the frequency pairs to specific windows
    modIdx(iWin,:,:,:) = modIdxTemp;

    % send(D, iWin); % Uncomment this line if you are tracking progress 
end

% Uncomment the following if you are tracking progress
% totalTime = toc;
% fprintf('Total time: %.1f minutes (%.2f hours)\n', totalTime/60, totalTime/3600); 

end

