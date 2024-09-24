function freq_centers = FFT_FreqBinCenters(NFFT,Fs) 
% Return the bin centers (Hz) for a given FFT, UDSP,pp 52 
nSpans = NFFT-1;
df     = Fs / NFFT;
freq_centers = df * (0:nSpans);
end