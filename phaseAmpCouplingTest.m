%% Simultate a few sinusoids


highFreqProbeA = pAGamma(1:1e3,:);
lowFreqWithin  = pATheta(1:1e3,:);
lowFreqBetween = pBTheta(1:1e3,:); 

fs = 1e3;
%%
for iType  = 1:4
     highFreq = pBGamma(1:1e3,:);

    switch iType
        case 1
            lowFreq  = lowFreqWithin;
            typeLabel = 'Within Probe';
        case 2
            lowFreq  = pATheta(1:1e3,:); 
            typeLabel = 'Between Probe - Theta';
        case 3
            lowFreq   = pAAlpha(1:1e3,:);
            typeLabel = 'Between Probe - Alpha';
        case 4
            lowFreq   = pABeta(1:1e3,:);
            typeLabel = 'Between Probe - Beta';
    end

    clear amplitude phase zComplex edges bin meanZ normMean ...
    zUni highFreqUni lowFreqUni
    minChannels = min([size(highFreq,2) size(lowFreq,2)]); 
    % Get analytical signal
    amplitude = abs(hilbert(highFreq(:,1:minChannels))); % Gamma
    phase     = angle(hilbert(lowFreq(:,1:minChannels))); % Theta

    % Complex valued signal
    zComplex = amplitude.*exp(1i*phase);

    [~,edges,bin] = histcounts(rad2deg(phase),'BinMethod','integers','BinEdges',-180:30:180);
    %
    for iMean = 1:size(zComplex,2)
        meanZ(:,iMean) = accumarray(bin(:,iMean),zComplex(:,iMean),[],@mean);
        % normMean = meanZ-mean(meanZ)./std(meanZ);
        normMean(:,iMean)  = abs(meanZ(:,iMean))./sum(abs(meanZ(:,iMean)));
    end

    % Plot the phase and amplitude in a polar plot
    figure; subplot(121); polarplot(median(phase,2,'omitnan'),median(amplitude,2,'omitnan'),'.');

    % Plot the distribution of means as a function of phase angle
    subplot(122); b =bar(edges(1:end-1)+180,median(normMean,2,'omitnan'),'histc');ylim([0 0.1]);
    b.FaceColor = [0 0.4470 0.7410]; 
    xlim([-20 380]); xlabel('Phase'); ylabel('Mean PAC (normalized)');box off;
    sgtitle(typeLabel);


    % Uniform distribution
    dt = 1/fs;
    t = (0:dt:1-dt)';
    highFreqUni = mean(abs(normMean)).*sin(2*pi*45*t);
    lowFreqUni  = mean(abs(normMean)).*sin(2*pi*5*t);

    ampUni   = abs(hilbert(highFreqUni));
    phaseUni = angle(hilbert(lowFreqUni));

    zUni = ampUni.*exp(1i*phase);
    meanZNew = accumarray(bin(:),zUni(:),[],@mean);
    normMeanUni = abs(meanZNew)./sum(abs(meanZNew));% meanZNew-mean(meanZNew)./std(meanZNew);

    figure;
    subplot(121); polarplot(phaseUni,ampUni,'.');

    subplot(122); b =bar(edges(1:end-1)+180,abs(normMeanUni),'histc');
    b.FaceColor = [0 0.4470 0.7410]; ylim([0 0.1]);
    xlim([-20 380]); xlabel('Phase'); ylabel('Mean PAC (normalized)');box off;
    sgtitle([typeLabel ' Uniform distribution']);

    % KL divergence
    modulationIdx(1,1:minChannels)       = sum(abs(normMeanUni).*log2(abs(normMeanUni)./abs(normMeanUni)),'omitnan')/log2(size(normMeanUni,1));
    modulationIdx(iType+1,1:minChannels) = sum(abs(normMean).*log2(abs(normMean)./abs(normMeanUni)),'omitnan')/log2(size(normMeanUni,1));

end
%%
figure; bar(abs(median(modulationIdx,2))); xticklabels({'Uniform','Within probe','theta-gamma','alpha-gamma','beta-gamma'}); box off;
ylabel('Modulation Index')

% Get the spectrogram? 
