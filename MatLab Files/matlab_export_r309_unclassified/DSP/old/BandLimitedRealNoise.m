function out = BandLimitedRealNoise(UniformNoise,Fs,Fc,BW)
    if ( (Fc + BW/2) > (Fs/2) )
        error('increase sample rate');
    end
    % first compute FFT
    fft_input = fft(UniformNoise);
    save_vector = zeros(1,length(fft_input));
    df = Fs / (length(fft_input)-1);
    fvec = 0:df:Fs;
    istart = find(fvec >= (Fc-BW/2),1,'first');
    istop  = find(fvec <= (Fc+BW/2),1,'last');
    save_vector(istart:istop) = 1;
    save_vector = save_vector(end:-1:1);
    save_vector(istart:istop) = 1;
    save_vector = save_vector(end:-1:1);
    
    Ascale = sum(save_vector)./length(save_vector);
    
    % to maintain same average power, divide by sqrt(Ascale)
    out = real(ifft(save_vector .* fft_input))./sqrt(Ascale);
end