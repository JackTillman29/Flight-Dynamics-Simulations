close all; clear all; clc

Fs = 500e6; dt = 1/Fs;

RF1 = 10e6;
RF2 = 10e6;
P_carrier_W = 1;
dur = 10e-6;

% dBc == decibels relative to carrier power
phs_noise_psd_dBc_Hz = -55 % dBc / Hz

t = [0:dt:dur-dt];
N = length(t);

%==========================================================================
% PHASE NOISE
%==========================================================================
phs_noise_pwr = P_carrier_W * 10^(phs_noise_psd_dBc_Hz / 10);
% phs_noise = sqrt(phs_noise_pwr) * 0.5 * (randn(1,N) + 1i*randn(1,N));
phs1_noise = sqrt(phs_noise_pwr) * randn(1,N);
phs2_noise = sqrt(phs_noise_pwr) * randn(1,N);

%==========================================================================
% NONLINEAR AMPLIFIER TRANSFER CURVE
%==========================================================================
asymm = 0.0;
amp_transfer_input = linspace(-100,100,10000);
amp_transfer_output = atan2d((amp_transfer_input+asymm),1) / 90;

figure; plot(amp_transfer_input,amp_transfer_output)

%==========================================================================
% DEFINE AMPLIFIER
%==========================================================================
preamp_gain1_dBV = -5;
preamp_gain2_dBV = -5;
preamp_gain1_V = 10^(preamp_gain1_dBV/10);
preamp_gain2_V = 10^(preamp_gain2_dBV/10);

%==========================================================================
% GENERATE SIGNALS
%==========================================================================
osc1_before_amp = real(exp(1i*2*pi*RF1.*t + 1i*phs1_noise));
osc1_after_amp  = interp1(amp_transfer_input,amp_transfer_output, preamp_gain1_V * osc1_before_amp);
osc2_before_amp = real(exp(1i*2*pi*RF2.*t + 1i*phs2_noise));
osc2_after_amp  = interp1(amp_transfer_input,amp_transfer_output, preamp_gain2_V * osc2_before_amp);

%==========================================================================
% MIX SIGNALS (AFTER SATURATION)
%==========================================================================
mixout_before_amp = osc1_before_amp .* osc2_before_amp;
mixout_after_amp = osc1_after_amp .* osc2_after_amp;

%==========================================================================
% COMPUTE SPECTRUMS
%==========================================================================
df = 1/dur;
frq = [-Fs/2:df:(Fs/2-df)];
fosc1_before_amp = fftshift(fft(osc1_before_amp));
fosc1_after_amp  = fftshift(fft(osc1_after_amp));
fosc2_before_amp = fftshift(fft(osc2_before_amp));
fosc2_after_amp  = fftshift(fft(osc2_after_amp));
fmixout_before_amp = fftshift(fft(mixout_before_amp));
fmixout_after_amp = fftshift(fft(mixout_after_amp));

%==========================================================================
% PLOT
%==========================================================================
if(~exist('hfig'))
    hfig = figure;
end
subplot(2,1,1);
plot(t*1e6,mixout_before_amp)
hold on;
plot(t*1e6,mixout_after_amp)
hold off;

colors = lines(2);
subplot(2,1,2);
hp(1) = plot(frq, 20*log10(abs(fmixout_before_amp)),'color',colors(2,:));
hold on;
hp(2) = plot(frq, 20*log10(abs(fmixout_after_amp)),'color',colors(1,:));
hold off;
hp = fliplr(hp);
set(gca,'children',hp)
legend(hp,'before','after')





