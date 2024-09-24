close all; clear all; clc

f1 = 1e6;
Fs = 20*f1; dt = 1/Fs;
t = [0:1000*(Fs/f1)-1]*dt;
t = [0:1000*(Fs/f1)]*dt;
s1 = exp(1i*2*pi*f1.*t);

f2 = 10e6;
s2 = exp(1i*2*pi*f2.*t);

s1 = s1 + s2;

fs1 = WaveObj(s1,[],Fs,length(s1),[],'z0',1);
plot(fs1)


fs2 = WaveObj(s1,[],Fs,length(s1),[],'z0',50);
% fs2 = SetZ0(fs2,50);
plot(fs2)

