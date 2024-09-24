close all;
clear all;
clc;

addpath('C:\Users\<user>\Documents\MATLAB\MATLAB\PulsesAndWaveforms');
addpath('C:\Users\<user>\Documents\MATLAB\MATLAB\Printing\');

t = 0:0.001:1;

y = SimpleLfm(25,125,0.5,0.001);
wo = WaveObj(y,0*y+1,1000,length(t));
plot(abs(fft(y)))

figure;
plot(wo);