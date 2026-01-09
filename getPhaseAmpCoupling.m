function modIdx = getPhaseAmpCoupling(phaseA,amplitudeA,phaseB,amplitudeB,winStart,winEnd,winSize,...
    chAIdx, chBIdx)
% Calculate the phase amplitude coupling.
% Based on Tort et al (2010), Journal of Neurophysiology 
% phaseVal: Low frequency phase
% amplitudeVal: High frequency amplitude
% Keerthana Manikandan
% October 13,2025
% Updated on Jan 9, 2026 

% Initialize variables
nLow     = size(phaseA,1);     % # Low frequency windows 
nHigh    = size(amplitudeA,1); % # High frequency windows
nWin     = numel(winStart);    % Windows to iterate
nChan    = size(chAIdx,1);     % Channel combinations

nBins    = 18;
binEdges = -180:20:180;
logN     = single(log(nBins)); 
groupIdx = single(repelem(1:nChan,winSize)');
phaseA   = rad2deg(phaseA);
pj       = zeros(nBins,nChan); %#ok<PREALL>

modIdx   = NaN(nWin,nHigh,nLow,nChan,'single');

% Setup progress tracking
D = parallel.pool.DataQueue;
progress = 0;
startTime = tic;

afterEach(D, @updateProgress);

    function updateProgress(~)
        progress = progress + 1;
        elapsed  = toc(startTime);
        pct = 100 * progress / nWin;
        avgTime = elapsed / progress;
        remaining = (nWin - progress) * avgTime;
        
        % Clear previous line and print update
        % fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        fprintf('Progress: %d/%d (%.1f%%) | Elapsed: %.1fm | Remaining: %.1fm | Avg: %.1fs/win', ...
            progress, nWin, pct, elapsed/60, remaining/60, avgTime);
    end

parfor iWin = 1:nWin 

    amplitudeATemp = squeeze(amplitudeA(:,winStart(iWin):winEnd(iWin),chAIdx)); %#ok<PFBNS>
    phaseATemp     = squeeze(phaseA(:,winStart(iWin):winEnd(iWin),chAIdx)); %#ok<PFBNS>

    amplitudeBTemp = squeeze(amplitudeB(:,winStart(iWin):winEnd(iWin),chBIdx)); %#ok<PFBNS>
    phaseBTemp     = squeeze(phaseB(:,winStart(iWin):winEnd(iWin),chBIdx)); %#ok<PFBNS>

    [~,~,binA]     = histcounts(phaseATemp,'BinEdges',binEdges); binA = single(binA);
    [~,~,binB]     = histcounts(phaseBTemp,'BinEdges',binEdges); binB = single(binB); 
    
    modIdxTemp    = NaN(nHigh, nLow, nChan, 'single');

    for iHigh = 1:nHigh
        for iLow = 1:nLow           

            ampAHigh = squeeze(amplitudeATemp(iHigh,:,:));
            binALow  = squeeze(binA(iLow,:,:));

            ampBHigh = squeeze(amplitudeBTemp(iHigh,:,:));
            binBLow   = squeeze(binB(iLow,:,:));

            % Reshape 'bin' and 'amplitude' into column vectors
            binCol       = binALow(:);
            amplitudeCol = ampAHigh(:);

            % % Divide phase into bins
            % [~,edges,bin] = histcounts(rad2deg(phaseVal),'BinEdges',-180:20:180);
            %
            % % Get the mean powers for each bin
            % n = size(edges,2)-1;
            %
            % % Reshape 'bin' and 'amplitude' into column vectors
            % binCol       = bin(:);
            % amplitudeCol = amplitudeVal(:);
            %
            % % Create a grouping index for each column
            % numColumns = size(bin, 2);
            % numRows = size(bin, 1);

            pj = accumarray([binCol, groupIdx], amplitudeCol, [], @mean);

            % pj = mean(amplitude) at phase j

            pj = pj./sum(pj,1,'omitnan'); % normalized pac = pj/sum(pj for j= 1:n) % The sum should be equal to 1

            % Calculate modulation index
            modIdxTemp(iHigh, iLow, :)  = (logN-(-sum(pj.*log(pj))))./logN;
        end
    end
    modIdx(iWin,:,:,:) = modIdxTemp;
    send(D, iWin);
end

totalTime = toc;
fprintf('Total time: %.1f minutes (%.2f hours)\n', totalTime/60, totalTime/3600); 

end

