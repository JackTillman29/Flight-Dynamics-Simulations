function out = BandLimitedComplexNoise3(N,Fs,Fc,BW,PSD)
    if ( ( Fc + BW/2) > Fs )
        error('increase sample rate');
    end

    % if 1-dimensional
    if(any(size(N) == 1))
        N = [1 N];
    else
        M = N;
        N = 1;
    end
    UniformNoise = sqrt(PSD*Fs) * exp(1i*2*pi*rand(N));
    
    % First compute FFT
    fft_input = fft(UniformNoise);
    save_vector = zeros(1,length(fft_input));
    df = Fs / (length(fft_input));
    %disp(['DF / BIN = ' num2str(df) 'Hz']);
    dt = 1 / Fs;
    fvec = 0:df:Fs;
    if(BW ~= Fs)
        istart = find(fvec >= (Fc - BW/2),1,'first');
        istop = find(fvec <= (Fc+BW/2),1,'last');
    else
        istart = 1;
        istop = length(fvec)-1;
    end
    window = 0*hamming(length((istart:istop)))+1;
    save_vector(istart:istop) = window;
    
    %Ascale = sum(save_vector)./length(save_vector);
    
    % out = ifft(save_vector .* fft_input)./sqrt(Ascale);
    out = ifft(save_vector .* fft_input);
    
    if(nargout==0 & any(size(N) == 1))
        %t = 0:(length(fvec)-1);
        %t = t .* dt;
        x = WaveFftStruct(out,0*out+1,Fs);
        figure;
        subplot(3,1,1);ah = gca;
        plot(x.frq,10*log10(abs(x.fft)));grid on;
        subplot(3,1,2);ah = [ah gca];
        plot(x.frq,10*log10([ ...
            abs(x.psd).' ...
            ]));grid on;
        subplot(3,1,3);ah = [ah gca];
        plot(x.frq,(abs(x.cpsd)));grid on;
        linkaxes(ah,'x');

        % Total Power
        PT = BW * PSD;
        out2 = ifft(sqrt(PT * df / BW) * N * save_vector);
        x = WaveFftStruct(out2,0*out2+1,Fs);
        
        
        figure;
        subplot(3,1,1);ah = gca;
        plot(x.frq,10*log10(abs(x.fft)));grid on;
        subplot(3,1,2);ah = [ah gca];
        plot(x.frq,10*log10(abs(x.psd)));grid on;
        subplot(3,1,3);ah = [ah gca];
        plot(x.frq,(abs(x.cpsd)));grid on;
        linkaxes(ah,'x');
    end
    
end