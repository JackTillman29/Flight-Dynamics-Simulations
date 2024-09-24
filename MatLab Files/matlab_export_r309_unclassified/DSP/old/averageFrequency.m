function [sigAvgFreq] = averageFrequency(s,Fs)
% computes the average frequency and the signal BW of the input signal

%% average freq calculation
N = length(s);
dt = 1/Fs;
% ws = fftshift(fft(s,N)/N);
% df = Fs / N;
% f = [-Fs/2:df:Fs/2-df];
% dw = df*2*pi/Fs;
% w = f*2*pi/Fs; % [-pi,pi)
% 
% eds = abs(ws).^2; % energy density spectrum
% % psd = eds / df;
% % cpsd = cumsum(psd)*df;
% sigPower = sum(eds);
% 
% % should be equal to each other (Parseval's / Rayleigh's theorem)
% % energy_t = (1/N)*sum(abs(s).^2)
% % energy_f = sum(abs(ws).^2)
% 
% avgWfn = cumsum(w.*abs(ws).^2)/sigPower; % NEED TO SCALE BY TOTAL SIGNAL POWER
% avgFfn = avgWfn * Fs/(2*pi); % rescale back to frequency from radians
% sigAvgFreq = avgFfn(end);



% to calculate the average frequency of a signal:
%  integrate( conj(s) * (1/j) * d/dt(s) * dt )
%    i.e. compute derivative of s then integrate according to above eq
dsdt = diff(s)/dt;
d2sdt2 = diff(diff(s))/(dt.^2);
approx = conj(s(2:end)) .* (1/1i) .* dsdt;
% avgWint = cumsum(approx*dt);
% avgF = abs(avgFint(end))/N/sigPower
avgW = abs(sum(approx*dt))/(N*dt);
sigAvgFreq = avgW / (2*pi);
%  TODO: explain why we divide by N*dt and 2*pi
%      1) divide by N*dt (signal duration)
%            eval integral requires normalization by signal duration (N*dt)
%      2) divide by 2*pi
%            w = 2*pi*f, above eq is for average rad/sec

% sigBandwidth = 



%% bandwidth calculation



end