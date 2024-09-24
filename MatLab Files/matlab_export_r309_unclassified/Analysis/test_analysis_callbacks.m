close all; clear all; clc

addpath('I:\MATLAB\Signal Generation\LFM Signal\')

Fs = 10e6;
dt = 1/Fs;
N = 10000;
t = [0:(N-1)]*dt;

% figure
% add_analysis_callbacks
% x = randn(1,N);
% plot(t,x)

figure
add_analysis_callbacks;

IF = 1e6;
BW = 100e3;
PW = 250e-6;
t = [0:dt:PW];
N = length(t);
% x = exp(1i*2*pi*f.*t);
x = SimpleLfm(IF-BW/2,IF+BW/2,PW,dt);
plot(t,x)
