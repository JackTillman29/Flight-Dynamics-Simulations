function [mag,frq]=FFT_CmplxFilterBin(NFFT,Fs,WindowFcn,NoiseVoltage,TargetDF)
% Syntax: [mag,frq]=FFT_CmplxFilterBin(NFFT,Fs,WindowFcn,NoiseVoltage,TargetDF)
%         WindowFcn not necessary. Or it can be a fcn handle @hamming, etc
%         NoiseVoltage not necessary. Or it can be in Volts
%         if nargout = 0, automatically plots

    % need to sweep across a bin in frequency
    df = Fs/NFFT;
    dt = 1.0 / Fs;
    
    % granularity
    if(~exist('TargetDF'))
        sdf = df/100;
    else
        sdf = TargetDF;
    end
    
    % select a middle bin
    nBin = floor(NFFT/2);
    
    % sweep vector
    fvect = 0:sdf:Fs;
    tvect = (0:(NFFT-1)) * dt;
    
    % output magnitude
    nSamples = length(fvect);
    outputMag = zeros(3,nSamples);
    
    if ( nargin > 2 )
        if( isa(WindowFcn,'function_handle') )
            w = WindowFcn(NFFT);
        else
            disp('assuming direct window passed');
            w = WindowFcn;
        end
    else
        w = ones(1,NFFT);
        WindowFcn = @none;
    end
    
    if ( nargin > 3 )
        nz = NoiseVoltage * exp(1i*2*(rand(1,NFFT)-0.5)*pi);
    else
        nz = zeros(1,NFFT);
    end
    
    
    for k = 1:length(fvect) % for each frequency
        try
            signal = w .* exp(1i*2*pi*fvect(k).*tvect) + nz;
        catch
            signal = w' .* exp(1i*2*pi*fvect(k).*tvect) + nz;
        end
        signal_fft = fft(signal);
        outputMag(:,k) = signal_fft((nBin-1):(nBin+1));
    end % for each frequency
    
    % output mapping
    mag = outputMag ./ NFFT;
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