close all;
clear all;
clc;

f = 11700;

dt1 = 0.001;
dt2 = 10e-6;

t1 = 0:dt1:1;
t2 = 0:dt2:1;

y1 = exp(1i*2*pi*f*t1);
y2 = exp(1i*2*pi*f*t2);

figure;
plot(t1,real(y1),'.-',t2,real(y2),'o-');
legend('Low','High');

figure;
fft1 = WaveFftStruct(y1,0*t1+1,1000);
fft2 = WaveFftStruct(y2,0*t2+1,100000);

plot(fft1.frq,fft1.psd);
hold on;
plot(fft2.frq,fft2.psd);
hold off;