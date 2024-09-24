close all; clear all; clc

addpath('I:\MATLAB\Signal Generation\Trivial\')

Fs = 100e3; dt = 1/Fs;

T = 1;
t = [0:dt:T-dt];
x = sin(2*pi*100*t);

fclk = 367;
clk = impulseOsc(fclk,T,dt);

% figure; plot(clk)
% 
% return

y = sampleAndHold(x,clk);


figure;
plot(t,x,'LineWidth',2)
hold on;
plot(t,y)
xlim([0 0.05])

%% you can also pass in multi-channel signal matrices and a common clock
% have all channels be sampled-and-held

t = [0:dt:T-2*dt];

numChannels = 4;

% x = rand(numChannels,N);
f = repmat([20 100 200 1000].',[1 length(t)]);
t = repmat(t,[numChannels 1]);
x = sin(2*pi*f.*t);
clk = clockOsc(248,T,0.2,dt);
% clk = rand(1,length(t)) > 0.9;

y = sampleAndHold(x,clk);

figure;
% plot(t.',x.');
% hold on;
plot(t.',clk,'k');
hold on;
plot(t.',y.','--')
xlim([0 0.05])
