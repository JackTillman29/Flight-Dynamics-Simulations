function [out,rms] = BandLimitedComplexNoise2(N,Fs,Fc,NoiseBW,PSD,RcvrBW,WindowFcn)
  
    if ( ( Fc + NoiseBW/2) > Fs )
        error('increase sample rate');
    end
    
    UniformNoise = sqrt(PSD*Fs) * exp(1i*2*pi*rand(1,N));
    
    % First compute FFT
    fft_input = fft(UniformNoise);
    save_vector = zeros(1,length(fft_input));
    df         = Fs / (length(fft_input));
    fvec       = 0:df:Fs;
    istart     = find(fvec >= (Fc - RcvrBW/2),1,'first');
    istop      = find(fvec <= (Fc+RcvrBW/2),1,'last');
    
    
%   %apply windowing if requested   
    if ( nargin > 6 )
        w = WindowFcn(length(istart:istop));
    else
        w = 1;
        WindowFcn = @none;
    end
    save_vector(istart:istop) = w;
    
    %Ascale = sum(save_vector)./length(save_vector);
    
    % out = ifft(save_vector .* fft_input)./sqrt(Ascale);
%     PT = BW * PSD;
%     out = ifft(sqrt(PT * df / BW) * N * save_vector);
    out = ifft(save_vector .* fft_input);
    rms = sqrt(mean(abs(out).^2));
    
    if(nargout==0)
        %t = 0:(length(fvec)-1);
        %t = t .* dt;
        x = WaveFftStruct(out,0*out+1,Fs);
        figure;
        subplot(3,1,1);ah = gca;
        plot(x.frq,10*log10(abs(x.fft)));grid on;
        title('FFT Output');
        xlabel('Frequency (Hz)');ylabel('dBV');
        subplot(3,1,2);ah = [ah gca];
        plot(x.frq,10*log10(abs(x.psd)));grid on;
        title('PSD Output');
        xlabel('Frequency (Hz)');ylabel('dBW/Hz');
        subplot(3,1,3);ah = [ah gca];
        plot(x.frq,(abs(x.cpsd)));grid on;
        xlabel('Frequency (Hz)');ylabel('W');
        title('Cumulative Power Distribution');
        linkaxes(ah,'x');
    end
    
    
end