close all; clear all; clc

addpath('I:\MATLAB\Analysis\')
addpath('I:\MATLAB\Printing\')
addpath('I:\MATLAB\DSP\')
addpath('I:\MATLAB\Utilities\')
addpath('I:\MATLAB\Windows\')
addpath('I:\MATLAB\Signal Generation\LFM Signal\')

Fs = 60e6;
dt = 1/Fs;

IF = 20e6;
BW = 10e6;
BW = 0e6;
PW = 1e-3;

x = SimpleLfm(IF-BW/2,IF+BW/2,PW,dt);
% x2 = SimpleLfm(IF+0.5*BW/2,IF-0.5*BW/2,PW,dt);
x2 = 70*randn(1,length(x));
x = x + x2;
% NN = 100e3;
% x = [x zeros(1,NN) x zeros(1,NN) x zeros(1,NN)];
figure; plot(real(x))

% return
tmp = WaveObj(x,x*0+1,Fs,length(x));

plot(tmp)
L = 1000;
overlap = 0;

% % no window
% stft(tmp,L,[],overlap,L)
% caxis([-50 0])
% 
% % hamming window
stft(tmp,L,@hamming,overlap,L)
% caxis([-50 0])

% % % hamming window with strings and scale factors added
% % timefac = 1e3;
% % timestr = 'ms';
% % freqfac = 1e-6;
% % freqstr = 'MHz';
% % stft(tmp,L,@hamming,overlap,L,...
% %     timefac,timestr,freqfac,freqstr)
% % caxis([-50 0])


taft(tmp,L,@hamming,overlap,L,1e-6,'MHz')
