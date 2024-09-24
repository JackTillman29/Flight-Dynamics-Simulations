close all;
clear;
clc;

DS1054Z = instrfind('ModelCode','0x04CE');
if(isempty(DS1054Z))
DS1054Z = visa( ...
    'ni', ...
    'USB0::0x1AB1::0x04CE::DS1ZA163655269::INSTR');
end


%%


data = DS1054Z_FullCapture(DS1054Z,24e6);

figure;
plot(DS1054Z_Convert(data.ch1));

