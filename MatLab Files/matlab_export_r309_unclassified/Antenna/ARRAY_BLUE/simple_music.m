%simple_music
close all;clear all;clc;
F = 5e6;
wavelen = 3e8/F;

addpath('I:\MATLAB\Antenna\ARRAY_BLUE');

seasonde = array('SeaSonde');
seasonde.CreateRectangularArray(wavelen,1.0,.001,1,2,1);%0.001
seasonde.upd_redefineCoordinates();
seasonde.upd_elementGaindBV = 1.0;

Fs = 1e6;
Ts = 1.0 / Fs;

t = single(0:Ts:(1000*Ts));

% generate simple noise signal (just to verify angle determination)
signal = randn(1,length(t)) + 1i*randn(1,length(t));
%signal = exp(1i*2*pi*2000*t);

fromAz = 5*pi/180;
fromEl = 0;

array_signals = seasonde.upd_genIQ_Channels(signal,fromAz,fromEl);

% figure;
% plot(real(array_signals.'))
% seasonde.upd_setSteer(fromAz,0);
% hold on;
% plot(real(exp(1i*seasonde.phs_steer).' * array_signals).')
% 


% Form correlation matrix [2 x 2]
Rest = array_signals*array_signals'./length(t);

% 
[Vr,Ur]=eig(Rest);

% If there is only one dominant eigenvalue, the associated eigenvector is
% the steering vector towards that source
[~,idxLargest] = sort(diag(Ur),'descend');



% Compute the steering vector required to generate this eigenvector
seasonde.upd_getSteer(-angle(Vr(:,idxLargest(1)))) * 180/pi






