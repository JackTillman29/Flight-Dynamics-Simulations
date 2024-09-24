function [fftsig,freq,sumsig,tdrms,fdrms] = FFT_RealInput(TimeDomainWaveform,Fs,NFFT)
% [fftsig,freq,sumsig,tdrms,fdrms] = FFT_RealInput(TimeDomainWaveform,Fs,NFFT(opt))
%
% fftsig is in units of input units
% sumsig is in power space (input units^2)
% abs(fftsig.^2) = peak power (peak voltage squared)
% 0.5 * abs(fftsig.^2) = average power (rms voltage squared)
%
% K Sawmiller, Booz | Allen | Hamilton

    if(~exist('NFFT'))
        NFFT = length(TimeDomainWaveform);
    end
    
    % Compute frequency bin width
    df = Fs / (NFFT-1);
    
    % Compute Bin Center Vector (Hz)
    freq = 0:df:Fs;
    
    % Compute the time domain rms signal level (input space)
    tdrms  = sqrt(mean(TimeDomainWaveform.^2));
    
    % Compute the Fourier Transform (input space)
    % NOTE: Scaling must be performed if zero padding -> sqrt() term
    fftsig = fft(TimeDomainWaveform,NFFT) ./ sqrt(NFFT * length(TimeDomainWaveform));
    
    % Compute the average power running sum (power space)
    sumsig = cumsum(abs(fftsig.^2));
    
    % half fft used for real data (other half wrapped mirror), hence:
    HFFT = floor(NFFT/2);
    fftsig = 2.0 * fftsig(1:HFFT);
    sumsig = 2.0 * sumsig(1:HFFT);
    freq   = freq(1:HFFT);
    
    % error checking
    fdrms = sqrt(sumsig(end));
    %disp(['RMS Signal Comparison within ' num2str(100*(abs(tdrms-fdrms)/(tdrms))) '%']);
    
    %disp(['Spectral frequency resolution = ' num2str(df)]);
    
end