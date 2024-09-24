function [mag,frq]=FFT_SingleBinOutput(NFFT,Fs,WindowFcn,OutputBin,OutputDf)
% Syntax: [mag,frq]=FFT_CmplxFilterBin(NFFT,Fs,WindowFcn,NoiseVoltage,TargetDF)
%         WindowFcn not necessary. Or it can be a fcn handle @hamming, etc
%         NoiseVoltage not necessary. Or it can be in Volts
%         if nargout = 0, automatically plots

    % need to sweep across a bin in frequency
    df = Fs/NFFT;
    dt = 1.0 / Fs;
    
    % sweep vector
    fvect = 0:OutputDf:Fs;
    nvect = (0:(NFFT-1));
    tvect = nvect * dt;
    
    % output magnitude
    nSamples = length(fvect);
    outputMag = zeros(1,nSamples);
    
    % address window functions
    if ( nargin > 2 )
        w = WindowFcn(NFFT);
    else
        w = ones(1,NFFT);
        WindowFcn = @none;
    end
    
    % DFT bin "m" response is equal to (1/NFFT) * SUM
    % (x(n)*e^(-j2pi*m*n/N)) from 0 to N-1
    
    twopij = 2*pi*1j;
    m = OutputBin;
    InvNFFT = 1.0 / NFFT;
    
    for k = 1:length(fvect) % for each frequency
        try
            %signal = w .* exp(1i*2*pi*fvect(k).*tvect);
            signal = w .* exp( twopij .* (fvect(k) .* tvect) - m .* nvect ./ NFFT);
        catch
            %signal = w' .* exp(1i*2*pi*fvect(k).*tvect);
            signal = w' .* exp( twopij .* (fvect(k) .* tvect) - m .* nvect ./ NFFT);
        end
        %signal_fft = fft(signal);
        %outputMag(:,k) = signal_fft((nBin-1):(nBin+1));
        outputMag(k) = InvNFFT .* sum(signal);
        
    end % for each frequency
    
    % output mapping
    mag = outputMag ./ NFFT;
    nBin=OutputBin;
    frq = fvect - (nBin-1)*df;
    
    if ( nargout == 0 )
        figure;
        plot(frq*1e-3,20*log10(abs(mag)));
        grid on;
        xlabel('Offset Frequency (kHz)');
        ylabel('Power Attenuation (dB)');
        title({ ...
            ['FFT Filter Representation NFFT=' num2str(NFFT) ', Fs = ' num2str(1e-3*Fs) ' kHz'], ...
            ['Window Function = ' func2str(WindowFcn)]});
        legend( ...
            ['Bin ' num2str(nBin-1)], ...
            ['Bin ' num2str(nBin  )], ...
            ['Bin ' num2str(nBin+1)] );
    end
    
end