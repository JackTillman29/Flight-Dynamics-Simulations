close all; clear all; clc

rootdir = 'I:\MATLAB\';
addpath([rootdir 'DSP'])

Fs = 10e6; dt = 1/Fs;
PW = 1;%74e-6;
t = 0:dt:PW-dt;
IF = 5e6;

fmSine = 100e3;
fmBW = 20e3;
% freqArray should be setup to be the length of the desired PW
freqArray = IF + fmBW/2*sin(2*pi*fmSine.*t);
% freqArray = IF + fmBW*sin(2*pi*fmSine.*t);
% freqArray = IF + fmBW*sin(2*pi*fmSine.*t) + fmBW*sin(2*pi*1.524*fmSine.*t + 2*pi*rand);

x = ArbFMOsc(freqArray,dt);

% Perform STFT on signal to verify nonlinear FM
wx = WaveObj(x,[],Fs,length(x));

nwin = 100;
winLen = floor(length(x)/nwin);
overlap = 0.5;
tSc = 1e6; tU = '[\mus]';
fSc = 1e-6; fU = '[MHz]';
NFFT = 10*winLen
stft(wx,winLen,[],overlap,NFFT,tSc,tU,fSc,fU)
xlim([0 30])
% overlay the frequency array
hold on
hp=plot(freqArray*fSc,t*tSc,'r')
legend(hp,'Desired FM')
