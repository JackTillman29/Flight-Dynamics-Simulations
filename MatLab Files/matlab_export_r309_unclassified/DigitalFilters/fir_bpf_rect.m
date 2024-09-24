function h_BPF = fir_bpf_rect(Fcenter,BW,Fs,N)
% compute filter taps of a rectangular bandpass filter approximation
% Syntax:
%   h_BPF = fir_bpf_rect(Fcenter,BW,Fs,N)
%       Fcenter: center frequency of band (Hz)
%       BW:      bandwidth (Hz)
%       Fs:      sampling frequency (Hz)
%       N:       number of taps
%       h_BPF:   output array of taps
%
% author: Jeff Hole (Booz Allen Hamilton) 2019-12

% compute the lowpass first, then mix it up to the bandpass
%   LPF prototype Fc should be BW/2

h_LPF = fir_lpf_rect(BW/2, Fs, N);

% compute n array
if(mod(N,2) == 1)
    n = [-(N-1)/2:(N-1)/2];
else
    n = [-(N/2):(N/2-1)];
end

fc = Fcenter / (Fs/2);
wc = pi * fc;
h_BPF = cos( wc * n ) .* h_LPF;



end