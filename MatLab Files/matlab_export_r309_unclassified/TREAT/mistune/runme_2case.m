close all;
clear;
clc;
matlab_root = 'C:\Users\<user>\Documents\MATLAB\matlab2018\';
if(~exist(matlab_root))
    matlab_root = 'C:\Users\<user>\Documents\MATLAB\matlab2018\';
end
addpath([matlab_root 'Signal Generation\Impulse Methods']);
addpath([matlab_root 'DSP']);
addpath([matlab_root 'Signal Generation\Noise Signal']);
addpath([matlab_root 'Signal Generation\PWM Signal']);
addpath([matlab_root 'Analysis']);
addpath([matlab_root 'DigitalFilters\Iowa Hills MATLAB']);
addpath([matlab_root 'Printing']);

pp = PrepForPrint();

% Create the input pulse train
Fs = 162.5e6;
Ts = 1/Fs;
CPI = 4e-4;
Fif = 30e6;
Ferr = 100000;
phs_offset = pi/4;

rf_in = PWM_Signal(1e-6,10e-6,Ts,floor(Fs*CPI),1);
t = (0:length(rf_in)-1) * Ts;
rf_in = exp(1i*2*pi*Fif*t) .* rf_in;
rf_in = rf_in + 0.1*exp(1i*2*pi*rand(1,length(t)));
pulse_present = abs(rf_in) > 0.1;


tg_out1 = exp(1i*(2*pi*(Fif+2000)*t+phs_offset));
tg_out2 = exp(1i*(2*pi*(Fif+-98000)*t+phs_offset));

% simulate scoring
sc1 = rf_in ./ tg_out1;
sc2 = rf_in ./ tg_out2;
sc1(pulse_present == 0) = nan;
sc2(pulse_present == 0) = nan;
