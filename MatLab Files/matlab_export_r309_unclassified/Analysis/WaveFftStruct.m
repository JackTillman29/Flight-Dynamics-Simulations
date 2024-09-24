function outputStructure = WaveFftStruct(tdSignal,tdWindow,Fs,NFFT,Option)
% function [sig_fft, sig_fft_freq, sig_psd, sig_cpsd] = WaveFft(tdSignal,tdWindow,Fs,['onesided'])
%
% Properties:
% For Real Signals
%
%    cumulative PSD should equal the squared rms value of the input
%    (average power)
%
% For Complex Signals
%    cumulative PSD should equal the mean squared value of the input
%
if(~exist('NFFT','var'))
    NFFT = length(tdSignal);
end

if(nargin == 5)
    [sig_fft, sig_fft_freq, sig_psd, sig_cpsd] = WaveFft(tdSignal,tdWindow,Fs,NFFT,Option);
else
    [sig_fft, sig_fft_freq, sig_psd, sig_cpsd] = WaveFft(tdSignal,tdWindow,Fs,NFFT);
end


if(nargout == 1)
    outputStructure.fft  = sig_fft;
    outputStructure.frq  = sig_fft_freq;
    outputStructure.psd  = sig_psd;
    outputStructure.cpsd = sig_cpsd;
else
    figure
    hax(1)=subplot(3,1,1)
    plot(sig_fft_freq,10*log10(abs(sig_fft)))
    ylabel('FFT')
    hax(2)=subplot(3,1,2)
    plot(sig_fft_freq,10*log10(sig_psd))
    ylabel('PSD')
    hax(3)=subplot(3,1,3)
    plot(sig_fft_freq,sig_cpsd)
    xlabel('Freq [Hz]')
    ylabel('Cum PSD')
    
    linkaxes(hax,'x')
end

end