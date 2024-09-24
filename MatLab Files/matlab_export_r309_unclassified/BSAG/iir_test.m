%% Test for IIR implementation
close all;
clear all;
clc;
addpath('z:\sawmillkd\MATLAB\DigitalFilters');

% Define a discrete space
dt = single(0.0001);
Fs = 1/dt;
t  = 0:dt:2;
y  = sin(2*pi*100*t);
y2 = sin(2*pi*100*t + pi);

% First create a filter
[filt.num,filt.den] = ChebyshevFilter( ...
    20/Fs, ... % Break
    0, ...      % 0 = Low Pass, 1 = High Pass
    0, ...      % Percent Ripple (0 to 29)
    2);         % Number of poles

dbode(filt.num,filt.den,1,dt,Fs/1000);
tic;
[yf] = applyIIR(filt.num,filt.den,y,'double');
method1 = toc;
fprintf('Method 1 took %8.3f seconds\n',method1);

figure;
plot(t,[y' yf']);
