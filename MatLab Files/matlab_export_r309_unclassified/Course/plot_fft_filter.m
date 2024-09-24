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


Fs = 1000;
Ts = 1./Fs;
t = 0:Ts:1;

frqa = linspace(0,Fs/2,10000);
maga = 0*frqa;
for frq = 1:length(frqa)
    y = cos(2*pi*frqa(frq)*t);
    fy = fft(y)./length(t);
    maga(frq) = 2*abs(fy(250));
end
figure;
plot(frqa,20*log10(maga));
xlabel('Frequency (Hz)');
ylabel('20log_{10}(A)');
add_print_callbacks;