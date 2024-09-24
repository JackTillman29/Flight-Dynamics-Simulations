function y = applyFilter(FLT,u)

% =========================================
% BRUTE FORCE METHOD (COMMENTED)
% =========================================

%     n     = length(u);
%     nTaps = length(tapGains);
%     
%     % first output is after nTaps input
%     nOutputs = n - nTaps + 1;
%     
%     % define the first set of input data that will fill the stages
%     firstInput = 1;
%     lastInput  = firstInput + nTaps - 1;
%     
%     y = zeros(size(u,1),size(u,2));
%     
%     for k = 1 : nOutputs
%         % generate output
%         output     = sum(tapGains .* u(lastInput:-1:firstInput));
%         
%          % assign output
%         y(nTaps-1+k) = output;
%         
%         % update the input data set
%         firstInput = firstInput + 1;
%         lastInput  = lastInput  + 1;
%         
%     end

% =========================================
% FAST CONVOLUTION METHOD
% =========================================
%y =conv(u,FLT.tapGains,'same');
nOut = length(u) + length(FLT.tapGains) - 1;
y = ifft( fft(u,nOut) .* fft(FLT.tapGains,nOut) );
nTooMany = nOut - length(u);

% KDS removed this shaving stuff because it causes time misalignment
%nShaveFront = ceil(nTooMany/2);
%nShaveEnd = nTooMany - nShaveFront;
y = y(1:(end-nTooMany));



if(nargout == 0)
    figure;
    subplot(2,1,1);
    plot([u' y'],'.-');
    grid on;
    axis tight;
    xlabel('Sample (#)');
    legend('Input','Output');
    subplot(2,1,2);
    [u_fft, u_fft_freq, u_psd, u_cpsd] = WaveFft(u,blackman(length(u)),FLT.Fs);
    [y_fft, y_fft_freq, y_psd, y_cpsd] = WaveFft(y,blackman(length(u)),FLT.Fs);
    plot(u_fft_freq,10*log10(abs([u_fft' y_fft'])),'.-');
    axis tight;
    xlabel('Frequency (Hz)');
    ylabel('10Log_{10}');
    grid on;
    legend('Input','Output');
end

end