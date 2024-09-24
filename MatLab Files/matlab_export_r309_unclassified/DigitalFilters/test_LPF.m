close all;
clear all;
clear classes;
clc;

addpath('Iowa Hills MATLAB');

Ts = 0.001;
t = 0:Ts:10;

ustep = 0*t;
ustep(100:end) = 1;
usin = sin(2*pi*1.0*t);

u = ustep .* usin;
dlmwrite('test_input.txt',u,'delimiter','\n');

u_simulink = [t' u'];

fil = biquad_stages('test_LPF_biquad_1000Hz_LP2Hz.txt');

y = fil.applyFilter(u);
sim('test_LPF_simulink');

p = load('out.txt');  % u fortran, y fortran, yc fortran


figure;
hold on;
plot(t,[u' y']);
plot(t, (y_simulink.signals.values),'r.');
plot(t, p(:,2),'bo');
plot(t, p(:,3),'bx');
xlabel('Time (sec)');
ylabel('Amplitude');
legend('Input','IH Matlab','Simulink','IH Fort','IHc Fort');