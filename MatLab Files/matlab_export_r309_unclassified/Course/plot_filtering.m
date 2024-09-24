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
addpath([matlab_root 'Analysis']);
addpath([matlab_root 'DigitalFilters\Iowa Hills MATLAB']);
addpath([matlab_root 'Printing']);
addpath([matlab_root 'Visualization']);
pp=PrepForPrint();

Fs = 10000;
Ts = 1/Fs;
t = 0:Ts:0.25;
y = cos(2*pi*25*t) + 0.5*cos(2*pi*1000*t) + 0.01*rand(1,length(t));

flt = biquad_stages('filter_demo.txt');
fltn = biquad_stages('filter_demo_notch.txt');
y2 = flt.applyFilter(y);
y3 = fltn.applyFilter(y);

wo = WaveObj(y,[],Fs);
plot(wo,1,'Hz',0,'pwr')
wo_t_f_plot(wo,gcf,'Original');
clear wo;

wo = WaveObj(y2,[],Fs);
plot(wo,1,'Hz',0,'pwr')
wo_t_f_plot(wo,gcf,'Low Pass');
clear wo;

wo = WaveObj(y3,[],Fs);
plot(wo,1,'Hz',0,'pwr')
wo_t_f_plot(wo,gcf,'Notch');
clear wo;

%%
plot(flt,10000,Fs,'Hz',[-2000 2000]);
PrepForPrint(gcf,pp);set(gcf,'Position',[869.0000   89.0000  453.2000  597.6000]);
kprint(gcf,'filter_demo_lpf.png',180);
plot(fltn,10000,Fs,'Hz',[-2000 2000]);
PrepForPrint(gcf,pp);set(gcf,'Position',[869.0000   89.0000  453.2000  597.6000]);
kprint(gcf,'filter_demo_notch.png',180);