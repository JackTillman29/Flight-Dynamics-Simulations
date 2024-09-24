function y = CentralConvolution(a,b,corrFlag)
    if ( nargin == 2 )
        corrFlag = 0;
    end

    a = a - mean(a);
    b = b - mean(b);
    % get FFT of both signals
    fa = fft(a);
    fb = fft(b);
    
    % multiply in frequency domain (convolution in time)
    if ( corrFlag == 1 )
        yt = ifft(fa.*conj(fb));
    else
        yt = ifft(fa.*fb);
    end
    
    %yt = real(xcf)/(sqrt(ACF1(1))*sqrt(ACF2(1)));

    
end