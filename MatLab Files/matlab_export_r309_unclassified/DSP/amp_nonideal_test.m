close all; clear all; clc

matlab_path = 'I:\MATLAB';
addpath([matlab_path '\DSP'])

ideal_amp = amp_nonideal;
myamp = amp_nonideal;
ra = runningAvg;

amp_gain_dB  = 30; % dB
amp_P1dB_dBm = -5; % dBm

myamp.setGain(amp_gain_dB,'dB')
myamp.setP1dB(amp_P1dB_dBm,'dBm')
ideal_amp.setGain(amp_gain_dB,'dB');

IF = 1e6;
Fs = 51e6; dt = 1/Fs;
dwell = 0.1e-3;
t = [0:dt:(dwell-dt)];
N = length(t);

input_power_dBm = -32;
noise_pwr_dBm   = -90;

x = sqrt(10^(0.1*(input_power_dBm-30))) * sin(2*pi*IF.*t);
x = x + 10^((noise_pwr_dBm-30)/20) .* randn(1,N);

y = myamp.applyGain(x);
y_ideal = ideal_amp.applyGain(x);


fx = abs(fft(x));
fy = abs(fft(y));

% PLOT INPUT/OUTPUT TIME SERIES
figure; add_print_callbacks;
hax(1)=subplot(2,1,1);
plot(t,x)
ylabel('Input')

hax(2)=subplot(2,1,2);
colors = lines(1);
plot(t,y_ideal,'-','color',0.5*ones(1,3));
hold on;
plot(t,y,'color',colors)
linkaxes(hax,'x')
xlim([0 100e-6])
ylabel('Output')
xlabel('Time [s]')
% 
% figure;
% hax(1) = subplot(2,1,1);
% plot(t,10*log10(abs(x)))
% hax(2) = subplot(2,1,2);
% plot(t,10*log10(abs(y)))
% linkaxes(hax,'x')

% x = ra.getdata();

% PLOT INPUT/OUTPUT SPECTRA
figure; add_print_callbacks;
df = Fs/N;
f = [0:df:(Fs-df)];

hax(1) = subplot(2,1,1)
plot(f*1e-6,10*log10(fx))
ylabel('Input')

hax(2) = subplot(2,1,2)
plot(f*1e-6,10*log10(fy))
ylabel('Output')
linkaxes(hax,'xy')
xlim([0 Fs/2]*1e-6)
xlabel('Freq [MHz]')

