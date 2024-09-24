close all; clear all; clc;

pf = @(x) 20*log10(abs(x));

N = 10e3;
x = randn(1,N);

peakwid1 = 5;
peakwid2 = 5;
del2 = -4*pi;
gain1 = 5;
gain2 = 3;

th = pi*linspace(-3,3,N)*peakwid1;
env = gain1*sin(th)./th;
env = env.^2;

th = pi*linspace(-3,3,N)*peakwid2;

env2 = gain2*sin(th-del2)./(th-del2);
env2 = env2.^2;

x = x + env + env2;


mmx = movmean(abs(x),400);

%==========================================================================
% define cfar parameters
nTrainingCells = 300;
nGuardCells    = 600;
numStdDev      = 10^(12/20);

[y,thr,det,sidethr] = adaptiveCfarProcessor(mmx,numStdDev,nTrainingCells,nGuardCells);

sidethr{1}(1:nTrainingCells+nGuardCells) = 0;
sidethr{2}(1:nTrainingCells+nGuardCells) = 0;
sidethr{1}(end-(nTrainingCells+nGuardCells):end) = 0;
sidethr{2}(end-(nTrainingCells+nGuardCells):end) = 0;

dets = (mmx > sidethr{3})*1;
y = mmx*0;
y(dets == 1) = mmx(dets == 1);

figure('Position',[520 64 1023 734]); add_print_callbacks;
hax(1) = subplot(2,1,1);
plot(pf(mmx))
hold on;
plot(pf(thr))
plot(pf(y),'k.')
plot(pf(sidethr{3}))
legend('x','adap','cfar')

hax(2) = subplot(2,1,2);
plot(pf(mmx))
hold on;
plot(pf(sidethr{1}))
plot(pf(sidethr{2}))
plot(pf(sidethr{3}))

linkaxes(hax,'x')

