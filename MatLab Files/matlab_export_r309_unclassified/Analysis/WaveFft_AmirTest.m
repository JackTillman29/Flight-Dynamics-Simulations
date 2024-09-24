function [sig_fft, sig_fft_freq, sig_psd, sig_cpsd] = WaveFft_AmirTest(tdSignal,tdWindow,Fs,NFFT,Option,varargin)
% function [sig_fft, sig_fft_freq, sig_psd, sig_cpsd] = WaveFft(tdSignal,tdWindow,Fs,['onesided'])
%
% Properties:
% For Real Signals
%    may use "onesided" option for plotting the spectrum of REAL SIGNALS
%    cumulative PSD should equal the squared rms value of the input
%    (average power)
%
% For Complex Signals
%    don't use "onesided" option for COMPLEX SIGNALS
%    cumulative PSD should equal the squared peak value of the input
%

z0 = 1;
for k = 1:2:length(varargin)
    if(strcmpi(varargin{k},'z0'))
        z0 = varargin{k+1};
    end
end


if(~exist('NFFT','var'))
    NFFT = length(tdSignal);
end

cpsdScale = NFFT ./ length(tdSignal);

%N = length(tdSignal);
df = Fs / NFFT;

window_cum_sum = cumsum(tdWindow);
window_frac    = window_cum_sum(end) / NFFT;

p_wnd = mean(tdWindow.^2); % this is the power factor associated with the window function
%disp([window_frac p_wnd]);

sig_fft = fftshift(fft(tdWindow.*tdSignal,NFFT) ./ NFFT ./ window_frac);
sig_fft_freq = Fs * ((NFFT/2:(NFFT/2-1)) ./ NFFT);
%sig_psd = (abs(sig_fft) .^ 2) / df / p_wnd;
sig_psd = (abs(window_frac*sig_fft) .^ 2) / df / p_wnd / z0;
sig_cpsd = cumsum(sig_psd) .* df .* cpsdScale;

if(exist('Option','var'))
    if( strcmp(Option,'onesided'))
        if (isreal(tdSignal))
            idxkeep = (0:floor(NFFT/2))+1;
            sig_fft      = 2.0 * sig_fft(idxkeep);
            sig_fft_freq = sig_fft_freq(idxkeep);
            sig_psd      = 2.0 * sig_psd(idxkeep);
            sig_cpsd     = 2.0 * sig_cpsd(idxkeep);
        else
            disp('Your signal is not real. One-sided option ignored.')
        end
    end
end


end