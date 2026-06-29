function [f,P1] = getFFT(sig,Fs,PLOT)
    L = length(sig);
%     T = 1/Fs;
%     t = (0:L-1)*T;
    f = Fs*(0:round(L/2))/L;

    Y = fft(sig);
    P2 = abs(Y/L);
    P1 = P2(1:round(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    if PLOT
        figure;
        plot(f,(P1));
        xlabel('Freq (Hz)');
        ylabel('Power');
        title(['FFT of ''',inputname(1),'''']);
        grid on;
    end
end