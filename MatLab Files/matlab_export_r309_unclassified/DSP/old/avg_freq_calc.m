close all; clear all; clc

format short eng

% IF = [15.0e6 23.0e6 24.0e6 25e6].'
IF = 20e6*rand(6,1);
% IF = 20e6*ones(10,1);
nTones = length(IF);
trueAvgFreq = sum(IF)/nTones;
trueBandwidth = max(IF)-min(IF);

a = ones(1,nTones);
PW = 100e-6;

% IF = 25e6;
Fs = 200*max(IF); dt = 1/Fs;

t = [0:dt:PW-dt];
N = length(t);

s = zeros(1,N);
for k = 1:nTones
    x = a(k).*exp(1i*2*pi*IF(k).*t);
    s = s + x;
end

IF = 15e6;
BW = 10e6;
x = SimpleLfm(IF-BW/2,IF+BW/2,PW,dt);
trueAvgFreq = IF;
trueBandwidth = BW;

noise = 1.0*exp(1i*rand(1,length(x)));
E_x = sum(abs(x).^2)*dt
E_n = sum(abs(noise).^2)*dt

s = x + noise;
E_s = sum(abs(s).^2)*dt

% return
figure; plot(real(s))


averageFrequency(s,Fs)

% return




%% calculate average frequency using the spectrum of s
ws = fftshift(fft(s,N)/N);
df = Fs / N; f = [-Fs/2:df:Fs/2-df];
dw = df*2*pi/Fs; w = f*2*pi/Fs; % [-pi,pi)

eds = abs(ws).^2; % energy density spectrum
psd = eds / df;
ceds = cumsum(eds);
cpsd = cumsum(psd)*df;
sigPower = cpsd(end)
% figure; plot(f,ceds); title('CEDS') % same as CPSD?
% figure; plot(f,cpsd); title('CPSD') % same as CEDS?

% should be equal to each other (Parseval's / Rayleigh's theorem)
energy_t = (1/N)*sum(abs(s).^2)
energy_f = sum(abs(ws).^2)


figure; subplot(2,1,1); plot(f,abs(ws)); title('signal spectrum')
subplot(2,1,2); plot(f,cpsd); title('CPSD') % same as CEDS?

avgWfn = cumsum(w.*abs(ws).^2)/sigPower; % NEED TO SCALE BY TOTAL SIGNAL POWER
avgFfn = avgWfn * Fs/(2*pi); % rescale back to frequency from radians
avgW = avgWfn(end);
avgF = avgFfn(end);

disp('Tones present:')
IF
disp('True average freq:')
disp(trueAvgFreq)
disp('True signal bandwidth:')
disp(trueBandwidth)
disp('Estimated average freq:')
disp(avgF)

% calculate average squared frequency
avgWfn2 = cumsum(w.^2.*abs(ws).^2)/sigPower;
avgFfn2 = avgWfn2 * Fs/(2*pi);
avgF2 = avgFfn2(end);

% sigBandwidth = sqrt(abs(avgF2 - avgF.^2))
% disp(['True average freq: ',num2str(trueAvgFreq),' Hz'])
% disp(['Estimated average freq: ',num2str(avgF),' Hz'])

% figure; plot(f,avgFfn); title('Average Freq Integral Function')

% return

%% calculate average frequency using time-domain signal only
% to calculate the average frequency of a signal:
%  integrate( conj(s) * (1/j) * d/dt(s) * dt )
%    i.e. compute derivative of s then integrate according to above eq
ds = diff(s);
dsdt = diff(s)/dt;
d2sdt2 = diff(diff(s))/(dt.^2);

% maxds = max(abs(dsdt))
% maxds2 = max(abs(d2sdt2))

% figure; plot(t,real(s),'b'); hold on; plot(t(2:end),real(dsdt),'r')

approx = conj(s(2:end)) .* (1/1i) .* dsdt;
approx = conj(s(2:end)) .* (1/1i) .* ds;
avgWint = abs(cumsum(approx)) ./ ([2:N]*dt);
avgFint = avgWint / (2*pi);

avgW = abs(sum(approx))/(N*dt);
avgF = avgW / (2*pi);
%  TODO: explain why we divide by N*dt and 2*pi
%      1) divide by N*dt (signal duration)
%            eval integral requires normalization by signal duration (N*dt)
%      2) divide by 2*pi
%            w = 2*pi*f, above eq is for average rad/sec

figure; plot(t(2:end),avgFint); hold on; plot([t(2) t(end)],[IF-BW/2 IF+BW/2],'r')

disp('Estimated average freq from time series alone:')
disp(avgF)


% avgF2 = sum(abs(dsdt).^2)*dt/N/sigPower

% avgF = abs(sum(conj(s(2:end)) .* (1/1i) .* dsdt)) / ((N-1)*sigPower)
% avgF2 = sum(abs(d2sdt2).^2)*dt*dt/sigPower
% 
% sqrt(avgF2)

% avgF = cumsum(conj(s(1:end-1)) .* (1/1i) .* dsdt) * dt;

% figure;
% subplot(3,1,1); plot(t,real(s),'b');
% subplot(3,1,2); plot(t(2:end),real(dsdt),'r')
% subplot(3,1,3); plot(t(2:end),abs(avgF))

