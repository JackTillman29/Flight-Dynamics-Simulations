close all; clear all; clc

Fs = 60e3; dt = 1/Fs;

f = 100;
T = 0.1;

f = linspace(100,200,floor(T/dt));

dcy = 0.5;
% x = rectOsc(f,T,dcy,dt);
x = triOsc(f,T,dt);

figure; plot(x)

ylim(1.1*[-1 1])


