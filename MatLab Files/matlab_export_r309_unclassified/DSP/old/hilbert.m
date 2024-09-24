function y = hilbert(x)
% function y = hilbert(x)
    xfft = fft(x);
    nfft = length(x);
    xfft((floor((nfft-1)/2)):end) = 0; % get rid of negative frequencies
    xfft(1) = 0.5*xfft(1);
    y = 2*ifft(xfft);
end