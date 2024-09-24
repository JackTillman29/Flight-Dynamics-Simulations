close all; clear all; clc

Fs = 100e6; dt = 1/Fs;

RF = 1e6;
P_carrier_W = 1;
dur = 1e-3;

% dBc == decibels relative to carrier power
phs_noise_psd_dBc_Hz = -10 % dBc / Hz

t = [0:dt:dur-dt];
N = length(t);

phs_noise_pwr = P_carrier_W * 10^(phs_noise_psd_dBc_Hz / 10);
% phs_noise = sqrt(phs_noise_pwr) * 0.5 * (randn(1,N) + 1i*randn(1,N));
phs_noise = sqrt(phs_noise_pwr) * randn(1,N);


%===========================================
% BUILD OSCILLATOR
%===========================================

%     exp( 1i*2*pi*f*t + 1i*phs )
osc = exp(1i*2*pi*RF.*t + 1i*phs_noise);

osc_phs_noise_only = exp(1i*phs_noise);

% phs_noise_pwr_meas_W = mean(abs(osc_phs_noise_only).^2);
% disp(['phase noise power (avg): ',num2str(10*log10(phs_noise_pwr_meas_W)),' dBW'])



% figure; plot(real(osc))
wosc = WaveObj(osc,[],Fs,N);
plot(wosc)

wosc = WaveObj(osc_phs_noise_only,[],Fs,N);
plot(wosc)


%% saturation











