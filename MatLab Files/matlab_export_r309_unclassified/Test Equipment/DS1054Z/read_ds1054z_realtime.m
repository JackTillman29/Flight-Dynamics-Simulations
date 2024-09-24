close all;
clear;
clc;


DS1054Z = visa( ...
    'ni', ...
    'USB0::0x1AB1::0x04CE::DS1ZA163655269::INSTR');

%%


DS1054Z_RealTimeCh1(DS1054Z);

figure;
plot(DS1054Z_Convert(data));

