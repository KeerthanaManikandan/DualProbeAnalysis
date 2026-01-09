function [modIdxA2A,modIdxB2B,modIdxAamp2Bphase,modIdxAphase2Bamp] = getPhaseAmpCoupling(phaseA,...
    amplitudeA,phaseB,amplitudeB,winStart,winEnd,winSize,chAIdx, chBIdx)

% Calculate the phase amplitude coupling.
% Based on Tort et al (2010), Journal of Neurophysiology 
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
phaseB   = rad2deg(phaseB);
pjA2A    = zeros(nBins,nChan); %#ok<PREALL>

modIdxA2A   = NaN(nWin,nHigh,nLow,nChan,'single');
modIdxB2B   = NaN(nWin,nHigh,nLow,nChan,'single');

modIdxAamp2Bphase   = NaN(nWin,nHigh,nLow,nChan,'single');
modIdxAphase2Bamp   = NaN(nWin,nHigh,nLow,nChan,'single');

% Setup progress tracking
D = parallel.pool.DataQueue;
progress = 0;
startTime = tic;
nChars    = 1;
afterEach(D, @updateProgress);

    function updateProgress(~)
        progress = progress + 1;
        elapsed  = toc(startTime);
        pct = 100 * progress / nWin;
        avgTime = elapsed / progress;
        remaining = (nWin - progress) * avgTime;
        
        % Clear previous line and print update
        fprintf(repmat('\b', 1, nChars));
        nChars = fprintf('\rProgress: %d/%d (%.1f%%) | Elapsed: %.1fm | Remaining: %.1fm | Avg: %.1fs/win', ...
            progress, nWin, pct, elapsed/60, remaining/60, avgTime);
    end

parfor iWin = 1:nWin 

    amplitudeATemp = squeeze(amplitudeA(:,winStart(iWin):winEnd(iWin),chAIdx)); 
    phaseATemp     = squeeze(phaseA(:,winStart(iWin):winEnd(iWin),chAIdx)); %#ok<*PFBNS>

    amplitudeBTemp = squeeze(amplitudeB(:,winStart(iWin):winEnd(iWin),chBIdx)); 
    phaseBTemp     = squeeze(phaseB(:,winStart(iWin):winEnd(iWin),chBIdx)); 

    [~,~,binA]     = histcounts(phaseATemp,'BinEdges',binEdges); binA = single(binA);
    [~,~,binB]     = histcounts(phaseBTemp,'BinEdges',binEdges); binB = single(binB); 
    
    modIdxTempA2A    = NaN(nHigh, nLow, nChan, 'single');
    modIdxTempB2B    = NaN(nHigh, nLow, nChan, 'single');

    modIdxTempAamp2Bphase = NaN(nHigh, nLow, nChan, 'single');
    modIdxTempAphase2Bamp = NaN(nHigh, nLow, nChan, 'single');

    
    for iHigh = 1:nHigh
        % Get the amplitude for iHigh 
        ampAHigh      = squeeze(amplitudeATemp(iHigh,:,:));
        ampBHigh      = squeeze(amplitudeBTemp(iHigh,:,:));
        amplitudeACol = ampAHigh(:);
        amplitudeBCol = ampBHigh(:);

        for iLow = 1:nLow
            % Get the phase for iLow frequency pairs
            binALow = squeeze(binA(iLow,:,:));
            binBLow = squeeze(binB(iLow,:,:));
            binACol = binALow(:);         
            binBCol = binBLow(:);           

            % Calculate modulation index between phase and amplitude
            % A - A
            pjA2A = accumarray([binACol, groupIdx], amplitudeACol, [], @mean);      
            pjA2A = pjA2A./sum(pjA2A,1,'omitnan'); % The sum should be equal to 1
            modIdxTempA2A(iHigh, iLow, :)  = (logN-(-sum(pjA2A.*log(pjA2A))))./logN;

            % B - B
            pjB2B = accumarray([binBCol, groupIdx], amplitudeBCol, [], @mean);
            pjB2B = pjB2B./sum(pjB2B,1,'omitnan'); % The sum should be equal to 1
            modIdxTempB2B(iHigh, iLow, :)  = (logN-(-sum(pjB2B.*log(pjB2B))))./logN;

            % A amplitude to B phase
            pjA2B = accumarray([binBCol, groupIdx], amplitudeACol, [], @mean);
            pjA2B = pjA2B./sum(pjA2B,1,'omitnan'); % The sum should be equal to 1
            modIdxTempAamp2Bphase(iHigh, iLow, :)  = (logN-(-sum(pjA2B.*log(pjA2B))))./logN;

            % A phase to B amplitude
            pjB2A = accumarray([binACol, groupIdx], amplitudeBCol, [], @mean);
            pjB2A = pjB2A./sum(pjB2A,1,'omitnan'); % The sum should be equal to 1
            modIdxTempAphase2Bamp(iHigh, iLow, :)  = (logN-(-sum(pjB2A.*log(pjB2A))))./logN;
        end
    end

 
    % Assign the frequency pairs to specific windows
    modIdxA2A(iWin,:,:,:) = modIdxTempA2A;
    modIdxB2B(iWin,:,:,:) = modIdxTempB2B;

    modIdxAamp2Bphase(iWin,:,:,:) = modIdxTempAamp2Bphase;
    modIdxAphase2Bamp(iWin,:,:,:) = modIdxTempAphase2Bamp;

    send(D, iWin);
end

totalTime = toc;
fprintf('Total time: %.1f minutes (%.2f hours)\n', totalTime/60, totalTime/3600); 

end

