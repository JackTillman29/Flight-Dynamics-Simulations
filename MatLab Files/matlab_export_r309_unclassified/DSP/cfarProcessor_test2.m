close all; clear all; clc

nTrain = [0 100];
nGuard = [0 0];

N = 1000;
x = randn(1,N);
x(N/2:end) = 3;

x = abs(x);

cfarThresh = 10^(5/10);
[y,thr,cfardet] = cfarProcessor(x,cfarThresh,nTrain,nGuard);

figure;
plot(x);
hold on;
plot(thr);

%====================
nTrain = [100 0];
nGuard = [0 0];

[y,thr,cfardet] = cfarProcessor(x,cfarThresh,nTrain,nGuard);
% subplot(2,1,2)
% plot(x);
% hold on;
plot(thr);

legend('Signal Envelope','Right Side CFAR','Left Side CFAR')





