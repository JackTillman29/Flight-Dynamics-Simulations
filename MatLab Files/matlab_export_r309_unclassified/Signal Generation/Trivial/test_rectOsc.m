close all; clear all; clc

f = 100;
T = 1;

Fs = 50e3; dt = 1/Fs;

i=0
for dutyCycle = linspace(0.01,0.99,1000)
    i = i + 1;
    
    x = rectOsc(f,T,dutyCycle,dt);
    
    mx(i) = mean(x);
    
end

figure; plot(mx)





