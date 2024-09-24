function h_LPF = fir_lpf_rect(Fc,Fs,N)
% compute filter taps of a rectangular lowpass filter approximation
% Syntax:
%   h_LPF = fir_lpf_rect(Fc,Fs,N)
%       Fc:     cutoff frequency of filter (Hz)
%       Fs:     sampling frequency (Hz)
%       N:      number of taps
%       h_LPF:  output array of taps
%
% author: Jeff Hole (Booz Allen Hamilton) 2019-12

% computes the inverse fourier series (time domain) approximation of a 
% rectangle function in the frequency domain, which we can compute directly
% using the sinc() function
%   H(f) = rect(f) --> h(t) = sinc(t)
%   ** this is assuming the rectangular function is centered at 0 Hz and
%   which is a lowpass filter (LPF) prototype

wc = pi * (Fc / (Fs/2));

% compute n array
if(mod(N,2) == 1)
    n = [-(N-1)/2:(N-1)/2];
else
    n = [-(N/2):(N/2-1)];
end

% tap computation for lowpass FIR filter prototype
h_LPF = (1./(wc*n)) .* (sin(wc.*n));
h_LPF(n == 0) = 1;

% compute mean-squared value and normalize h_LPF by it
rms_val = mean(h_LPF);
h_LPF = h_LPF ./ rms_val;


%     h_zeros = roots(h_LPF);
%     
%     th = linspace(0,2*pi,1000);
%     figure;
%     plot(cos(th),sin(th),'k--');
%     hold on; axis equal;
%     plot(real(h_zeros),imag(h_zeros),'o');
%     title('Pole-Zero Diagram')
%     xlabel('Re{roots(filter taps)}')
%     ylabel('Im{roots(filter taps)}')


end


