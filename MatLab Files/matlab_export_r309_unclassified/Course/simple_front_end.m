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

PRF = 100e3;
PD = 1e-6;
CPI = 4e-3;
rxLO = 20e6; % signals will be generated at a 10MHz carrier
rxRF = 8e9; % signal RF (for proper phase alignment)
kbolt = 1.38e-23;
Tref = 290;
Nf = 6;

Fs = 100e6;
rxIF = 10e6;
Ts = 1/Fs;
NPD = round(PD/Ts);

% create reference rectangular pulse @ rxLO MHz
txWaveform = exp(1i*2*pi*rxLO*(0:NPD-1)*Ts);

nSamples = round(CPI/Ts);

rxImpulse = zeros(1,nSamples);
gateImpulse = zeros(1,nSamples);

% target data
tgtRanges = [65000 75200 500];
tgtAmplitudes = [1*1.0e-6 1*1.0e-5 1*1.0e-5];
tgtDopplers = [14000 17000 1*0];

gateSignal = 1; % which signal above is the gate centered

generateTargetSkin;
generateRangeGate;


% Add Noise
rxSamples = rxSamples + GenerateUnfilteredNoise(length(rxSamples),Fs,kbolt*Tref*10^(0.1*Nf),'gaussian_iq');
% Create time vector
rxTime = (0:length(rxSamples)-1) * Ts;

test_point(rxSamples,Fs,'RF Input');

% Load baseband filters
flt_pw_bb = biquad_stages('simple_front_end_bb_pw_10M.txt');
flt_prf_bb = biquad_stages('simple_front_end_lp_prf.txt');
flt_clut_notch1_bb = biquad_stages('simple_front_end_notch_0PRF.txt');
flt_clut_notch2_bb = biquad_stages('simple_front_end_notch_1PRF.txt');
flt_prf_notch_bb = flt_clut_notch1_bb + flt_clut_notch2_bb;
flt_dop = biquad_stages('simple_front_end_dop.txt');
flt_pw_if = flt_pw_bb.shiftStages(rxIF,Fs);

% Downconvert
rx_1 = ddc(rxSamples,-rxLO+rxIF,Ts);
test_point(rx_1,Fs,'IF Input');


rx_2 = flt_pw_if.applyFilter(rx_1);
test_point(rx_2,Fs,'IF Filtered Input');

rx_3 = rx_2 .* gateSamples;
test_point(rx_3,Fs,'Range Gated Input');


rx_4 = flt_prf_notch_bb.applyFilter(ddc(rx_3,-rxIF,Ts));
rx_4 = ddc(rx_4,rxIF,Ts);
test_point(rx_4,Fs,'Gated + Clutter Notch');

rx_5 = flt_prf_bb.applyFilter(ddc(rx_4,-rxIF,Ts));
test_point(ddc(rx_5,rxIF,Ts),Fs,'Gated + Clutter Notch + PRF');

linkaxes(findobj('Tag','testpoint1'),'x')
linkaxes(findobj('Tag','testpoint2'),'x')

% reduce the sample rate for numeric stability
decrate = 100;
Fs2 = Fs / decrate;
Ts2 = 1 / Fs2;

rx2_5 = rx_5(1:decrate:end);



dopFilterBankCenters = (-2000:500:2000) + 14e3;

for dc = 1 : length(dopFilterBankCenters)
    %flt_dop_now=flt_dop.shiftStages(rxIF + dopFilterBankCenters(dc),Fs);
    rx_dop_fil = flt_dop.applyFilter(ddc(rx2_5,-dopFilterBankCenters(dc),Ts2));
    rx_dop_fil = ddc(rx_dop_fil,dopFilterBankCenters(dc),Ts2);
    if(dc == 1)
        output = zeros(length(dopFilterBankCenters),length(rx_dop_fil));
    end
    output(dc,:) = rx_dop_fil;
    if(dc == 5)
        test_point(rx_dop_fil,Fs2,'Dop');
    end
end
figure;
imagesc(abs(output));
colorbar;

%% Plot Filter Bank
if(0)
for dc = 1 : length(dopFilterBankCenters)
    rx_dop_fil = flt_dop.shiftStages(dopFilterBankCenters(dc),Fs2);
    if(dc == 1)
        plot(rx_dop_fil,10000,Fs2/1000,'kHz',([-3000 3000]+14e3)./1000);
        h1 = subplot(2,1,1);
        h2 = subplot(2,1,2);
        co = get(gca,'ColorOrder');
        co = [co;co;co;co];
    else
        plot(rx_dop_fil,10000,Fs2/1000,'Hz',([-3000 3000]+14e3)./1000);
        set(get(subplot(2,1,1),'Children'),'Parent',h1,'Color',co(dc,:));
        set(get(subplot(2,1,2),'Children'),'Parent',h2,'Color',co(dc,:));
        close(gcf);
        
    end
end

PrepForPrint(gcf,pp);
end
%% plot traces
axHandles = subplots_full(figure,1,length(dopFilterBankCenters),0);
for k = 1 : length(axHandles)
    plot(axHandles(k),real(output(k,:)),'Color',0.8*[1 0 0]);
    ht=text(10,0,[num2str(1e-3*dopFilterBankCenters(k),'%8.3f') 'KHz']);
    set(ht,'Fontsize',8,'Fontweight','bold');
    set(ht,'Parent',axHandles(k),'BackgroundColor','w');
end
set(axHandles,'Visible','on','XTick',[],'YTick',[],'XColor',0.9*[1 1 1])
linkaxes(axHandles,'y')
PrepForPrint(gcf,pp);
set(gcf,'Position',[931.4000   65.0000  307.2000  686.0000]);

%plot([flt_dop14 flt_dop14L flt_dop14H],100000,Fs*1e-3,'kHz',[0 40000]*1e-3)


