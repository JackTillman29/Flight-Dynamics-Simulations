close all; clear all; clc

addpath('I:\MATLAB\DigitalFilters\Iowa Hills MATLAB\')

Fs = 120e6;
% filt_filename = '6 Pole Inv Cheby LPF.txt';
% f = biquad_stages(filt_filename);
f = biquad_stages(3);

a0 =  1.000000000000000000;
a1 =  -1.733014396422086540;
a2 =  0.752347935115272315;
b0 =  0.301903596237685923;
b1 =  -0.584473653782186187;
b2 =  0.301903596237685923;

f = f.add_stage( biquad(b0,b1,b2, a1, a2) );

% 2.46 MHz LPF Fc
a0 =  1.000000000000000000;
a1 =  -1.818386036050307240;
a2 =  0.835887334603790055;
b0 =  0.148487430047835073;
b1 =  -0.279473561542187332;
b2 =  0.148487430047835073;

f = f.add_stage( biquad(b0,b1,b2, a1, a2) );

a0 =  1.000000000000000000;
a1 =  -1.927692091097539250;
a2 =  0.944006011814306545;
b0 =  0.022075930046734207;
b1 =  -0.027837939376701211;
b2 =  0.022075930046734207;

f = f.add_stage( biquad(b0,b1,b2, a1, a2) );


% duplicate filter
fUpper = f;
fLower = f;

i = 0;
for Fc = linspace(5e6,50e6,20)
    i = i + 1;
% shift filters up and down by equal amounts
% Fc = 40e6;
fUpper = fUpper.shiftStages(Fc,Fs);
fLower = fLower.shiftStages(-Fc,Fs);

% fUpper.plotPoleZero
% fUpper.plot
% fLower.plot

% combine filters
f_BPF = fUpper;
f_BPF = f_BPF.add_stage( fLower.biquads(1) ); f_BPF.nstages = f_BPF.nstages + 1;
f_BPF = f_BPF.add_stage( fLower.biquads(2) ); f_BPF.nstages = f_BPF.nstages + 1;
f_BPF = f_BPF.add_stage( fLower.biquads(3) ); f_BPF.nstages = f_BPF.nstages + 1;

% f_BPF.plotPoleZero
f_BPF.plot
% pause
kprint(gcf,['C:\Users\holeja\Desktop\figs\',num2str(i),'.png'],180);
close all

fUpper = fUpper.unshiftStages;
fLower = fLower.unshiftStages;

end






