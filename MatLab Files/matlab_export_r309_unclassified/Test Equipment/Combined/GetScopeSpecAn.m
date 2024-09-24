% Open up sockets

% close all active connections
% h=instrfind;
% for k = 1 : length(h)
%     if( strcmpi(h(k).Status,'open') )
%         fclose(h(k));
%     end
% end
addpath('C:\Users\Keith\Documents\MATLAB\MATLAB\Printing')
% reset all connectsions
instrreset;


SSA = visa('ni','USB0::0xF4EC::0xEE38::XXXXXXXXXX::INSTR');
SDG805 = visa('ni','USB0::0xF4ED::0xEE3A::SDG00003140697::INSTR');
DS1054Z = visa('ni','USB0::0x1AB1::0x04CE::DS1ZA163655269::INSTR');

DS1054Z.InputBufferSize = 65535;
DS1054Z.OutputBufferSize = 65535;
SSA.InputBufferSize = 65535;
SSA.OutputBufferSize = 65535;

SDG805.InputBufferSize = 65535;
SDG805.OutputBufferSize = 65535;

%% Set Up Sig Gen
SDG805_Sine(SDG805,30e6,0,'dbm');
SDG805_OutputImpedance(SDG805,'50');

% make sure we are in 50ohm impedance mode
imp = check_scpi_value(SDG805,'C1:OUTP?','LOAD');
if( ~strcmpi(imp,'50') )
    error('Not in 50Ohm Mode!!!');
end
SDG805_OnOff(SDG805,'on');

%%
pp = PrepForPrint();
pp.titleFontSize = 14;
[ssa_data,ssa_frq]=SSA3021X_GetScreenData(SSA,'dummy');
[scope_data] = DS1054Z_GetScreenData(DS1054Z,'dummy');

% NOTE: Calibrated test setup to include cable losses and coupler
% this way both time & frequency domain will align in amplitude
sf = 632/81.6;

figure;
subplot(1,2,1);
plot(scope_data.time*1e6,sf*1e3*scope_data.ch1);
xlabel('Time (\mus)');
ylabel('Volts (mV)');
grid on;
title({'Time Domain'});

subplot(1,2,2);
plot(1e-6*ssa_frq,ssa_data);
xlabel('Frequency (MHz)');
ylabel('db rel 0dBm');
grid on;
title('Frequency Domain');

PrepForPrint(gcf,pp);
set(gcf,'Color','w','Position', ...
     [-1637.8 770.2  1278  430.4]);



kprint(get(gcf,'Number'),'Demo_Pulsed_Zoom.png',180);