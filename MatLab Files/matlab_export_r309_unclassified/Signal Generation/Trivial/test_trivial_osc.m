close all; clear all; clc


Fs = 20e3;
dt = 1/Fs;
% dt = 0.00001;

f = 100;
T = 1;
t = [0:dt:T-dt];

zoom = 0.05;

%% build simple sawtooth oscillator
x = sawOsc(f,T,dt);

figure; plot(t,x); xlim([0 zoom*T])

%% difference of two saws gives a square oscillator with adjustable duty
% cycle
dutyCycle = 0.75;
x1 = sawOsc(f,T,dt,1-dutyCycle);
x2 = sawOsc(f,T,dt);
x = x1 - x2;

figure; plot(t,x); xlim([0 zoom*T]); title(['avg: ',num2str(mean(x))])

%% integrating a square wave gives a triangle wave
% close all;
dutyCycle = 0.5; % TODO: figure out a way so that any duty cycle will give 
                 %       a DC-compensated triangle waveform. (DC-compensated
                 %       means that the average of one complete cycle = 0)
xsq = rectOsc(f,T,dutyCycle,dt);

x = cumsum(xsq);

samplesPerPeriod = (1/f) / dt;
nRiseSamples = samplesPerPeriod / 2;
nFallSamples = samplesPerPeriod / 2;
% halfSamPerPeriod = samplesPerPeriod / 2;
x = 2 * cumsum(xsq) / nRiseSamples - 1; % divide by f to normalize.

figure; plot(x); xlim([0 1e3])
% figure; plot(t,x); xlim([0 zoom*T]); title(['avg: ',num2str(mean(x))])


%% clock oscillator
% close all;
dutyCycle = 0.2;
x = clockOsc(f,T,dutyCycle,dt);

figure; plot(t,x); xlim([0 zoom*T]);





