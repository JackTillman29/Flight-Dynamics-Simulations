close all; clear all; clc

matlab_path = 'I:\MATLAB';
addpath([matlab_path '\DSP'])

myamp = amp_nonideal;
ra = runningAvg;

amp_gain_dB  = 30; % dB
amp_P1dB_dBm = -5; % dBm

myamp.setGain(amp_gain_dB,'dB')
myamp.setP1dB(amp_P1dB_dBm,'dBm')

IF = 1e6;
Fs = 51e6; dt = 1/Fs;
dwell = 0.1e-3;
t = [0:dt:(dwell-dt)];
N = length(t);

input_power_dBm = -32;

noise_pwr_dBm = -90;

input_power_dBm_array = [-36:0.05:-10];
ntrade = length(input_power_dBm_array);
fx = zeros(ntrade,N);
fy = zeros(ntrade,N);

itrade = 0;
for input_power_dBm = input_power_dBm_array
    itrade = itrade + 1
x = sqrt(10^(0.1*(input_power_dBm-30))) * sin(2*pi*IF.*t);
x = x + 10^((noise_pwr_dBm-30)/20) .* randn(1,N);

% rax = rax.newdata(x);

y = myamp.applyGain(x);

fx(itrade,:) = abs(fft(x));
fy(itrade,:) = abs(fft(y));


% [junk,junk2,fx_dets] = cfarProcessor(fx(itrade,:),10^(4/10),10,3);
% [junk,junk2,fy_dets] = cfarProcessor(fy(itrade,:),10^(4/10),10,3);

a = 1;

% % PLOT INPUT/OUTPUT TIME SERIES
% figure;
% 
% hax(1)=subplot(2,1,1);
% plot(t,x)
% 
% hax(2)=subplot(2,1,2);
% plot(t,y)
% linkaxes(hax,'x')
% xlim([0 100e-6])
% 
% figure;
% hax(1) = subplot(2,1,1);
% plot(t,10*log10(abs(x)))
% hax(2) = subplot(2,1,2);
% plot(t,10*log10(abs(y)))
% linkaxes(hax,'x')

% x = ra.getdata();

% % % % PLOT INPUT/OUTPUT SPECTRA
% % % if(~exist('hfig_spectra'))
% % %     hfig_spectra = figure;
% % % end
% % % hax(1) = subplot(2,1,1)
% % % plot(t,10*log10(fx(itrade,:)))
% % % ylabel('Input')
% % % 
% % % hax(2) = subplot(2,1,2)
% % % plot(t,10*log10(fy(itrade,:)))
% % % ylabel('Output')
% % % linkaxes(hax,'x')
% % % 
% % % drawnow;
end


figure
imagesc((fy(:,1:round(N/2))));
