close all;
clear all;
clc;

addpath('I:\MATLAB\Analysis');

% BEFORE USING gen_pd_impulse(), PLEASE TAKE A LOOK AT THE M-CODE AND
% UNDERSTAND THE CAVEATS. THIS ROUTINE SHOWS HOW TO USE THIS FUNCTION.
% OTHER SCRIPTS SHOW CASES WHERE THESE CAVEATS ARE IMPORTANT TO UNDERSTAND.

PRF = 100e3;
dutyCycle = 0.3;
CPI = 5e-3;
rxLO = 17.16e6; % signals will be generated at a 10MHz carrier
rxRF = 10e9; % signal RF (for proper phase alignment)


Fs = 200e6;
Ts = 1/Fs;

PW = dutyCycle / PRF;
NPW = round(PW/Ts);

% create reference rectangular pulse @ 10MHz
txWaveform = exp(1i*2*pi*rxLO*(0:NPW-1)*Ts);

nSamples = round(CPI/Ts);

rxImpulse = zeros(1,nSamples);

% target data
tgtRange = 65000;
tgtAmplitude = 0.5;
tgtDoppler = 25e3;

% generate impulse return for this target
rxImpulse = gen_pd_impulse(...
    rxImpulse,...
    tgtRange,...
    PRF,...
    Ts,...
    tgtAmplitude,...
    rxLO,...
    rxRF,...
    tgtDoppler);

rxSamples = conv(rxImpulse,txWaveform);
rxTime = (0:length(rxSamples)-1) * Ts;

wo = WaveObj(rxSamples,[],Fs);
wo_ref = WaveObj(exp(1i*2*pi*rxLO*rxTime),[],Fs);

figure;
plot(rxTime,[exp(1i*2*pi*rxLO*rxTime).' rxSamples.' ])
title('Time domain');

figure;
%plot(20*log10(abs(fft([exp(1i*2*pi*rxLO*rxTime).' rxSamples.']))));
plot([wo wo_ref],1e-6,'MHz');

    
