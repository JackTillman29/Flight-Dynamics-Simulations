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

for ferr2 = 97e3; %linspace(-100e3,100e3,101)

tg_out = exp(1i*(2*pi*(Fif+ferr2)*t+phs_offset));

% simulate scoring
score_data = rf_in ./ tg_out;
score_data(pulse_present == 0) = nan;

plot(1e6*t,angle(score_data),'.');
xlabel('Time (\musec)');
ylabel('\Delta\phi (rad)');
title(['Tune Error: ' num2str(ferr2) ' Hz']);
ylim([-pi pi]);
PrepForPrint(get(gcf,'Number'),pp);

kprint(get(gcf,'Number'),['phase_plot_err_' num2str(ferr2) 'hz.png'],180);
close(gcf);
end