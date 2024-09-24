close all; clear all; clc

addpath('I:\MATLAB\Analysis');


%% CAVEAT: IN GENERAL, PULSE TRAIN IS NOT COHERENT WITH A FREE-RUNNING COHO
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



%==========================================================================
%  IMPLEMENT A NON-MOVING TARGET AND TEST FOR COHERENCY OVER THE DWELL
%==========================================================================

% target data
tgtRange = 0;
tgtAmplitude = 1;
tgtDoppler = 0;

% generate impulse return for this target
rxImpulse = zeros(1,nSamples);
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

figure;
plot(rxTime,[exp(1i*2*pi*rxLO*rxTime).' rxSamples.' ])
title('Time domain (non-moving target)');

%==========================================================================
%  IMPLEMENT A DOPPLER SHIFT AND TEST FOR COHERENCY AT THE DOPPLER SHIFT
%==========================================================================
% target data
tgtRange = 0;
tgtAmplitude = 1;
tgtDoppler = 12.437e3;

% generate impulse return for this target
rxImpulse = zeros(1,nSamples);
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

gate = 1*(abs(rxSamples)>0);

figure;
plot(rxTime,[exp(1i*2*pi*(rxLO+tgtDoppler)*rxTime).' rxSamples.' ])
title('Time domain (doppler-shifted target)');

error = gate.*exp(1i*2*pi*(rxLO+tgtDoppler)*rxTime) - rxSamples;

figure('Position',[21 85 1004 681]);
plot(rxTime,error)
title({'DEMO OF CAVEAT #1 (in gen_pd_impulse.m)';...
    'These errors are due to the fact that each pulse has a constant phase shift';...
    'So, in reality, the Doppler shift is only capture on a pulse-to-pulse';...
    'basis, not over the pulse duration of a single pulse (each pulse is';...
    'really at the IF (or rxLO) only!'},'Interpreter','none');


