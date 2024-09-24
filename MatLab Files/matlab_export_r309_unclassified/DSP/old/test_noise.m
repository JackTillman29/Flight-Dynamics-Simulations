close all;clear;clc;
pp = PrepForPrint;
Fs = 1e6;
BW = 30e3;
Fc = 100e3;

k = 1.3805e-23;
Tref = 290;
Nf = 10^(4/10);

PSD = k*Tref*Nf; % Watts/Hz

PAVG = PSD * BW;

Dur =3e-3;
N = floor(Dur*Fs);
t = (0:(N-1))*(1/Fs);

out = BandLimitedComplexNoise3(N,Fs,Fc,BW,PSD);

voltage = sqrt(2)*real(out);

outf =  WaveFftStruct(voltage,0*out+1,Fs,'onesided');

figure;
plot(t,1e6*voltage);
xlabel('Time (s)');
ylabel('\muV');
title({'Time History',['Vrms = ' num2str(std(voltage))]});
add_print_callbacks;

figure;
hist(1e6*voltage,100);
stddev = std(1e6*voltage);
xlabel('\muV');
ylabel('Occurences');
title({'Voltage Distribution',['\sigma=' num2str(stddev) ' \muV']});
add_print_callbacks;

figure;
hist(voltage.^2,100);
stddev = std(1e6*voltage.^2);
xlabel('\muV');
ylabel('Occurences');
title({'Voltage Distribution',['\sigma=' num2str(stddev) ' \muV']});
add_print_callbacks;

figure;
plot(outf.frq,outf.cpsd);
xlabel('Frequency (Hz)');
ylabel('W');
title('Cumulative Power');
add_print_callbacks;


%% loop for more data
nloops = 10000;
odat = zeros(nloops,2);
for jj = 1 : nloops
    out = BandLimitedComplexNoise3(N,Fs,Fc,BW,PSD);
    voltage = sqrt(2)*real(out);
    Vrms = std(voltage);
    Pavg = Vrms.^2;
    odat(jj,:) = [Vrms Pavg];
end
%%

mean_vrms = mean(odat(:,1));
mean_pavg = mean(odat(:,2));
std_vrms = std(odat(:,1));
std_pavg = std(odat(:,2));

figure;
subplot(2,2,1);
plot(1e9*odat(:,1));xlabel('CPI #');ylabel('Vrms (nV)');
title(['Vrms over ' num2str(nloops) ' CPIs']);

subplot(2,2,2);
plot(odat(:,2));xlabel('CPI #');ylabel('Pavg (W)');
title(['Pavg over ' num2str(nloops) ' CPIs']);

subplot(2,2,3);
hist(1e9*odat(:,1),50);xlabel('Vrms nV');ylabel('# Occurences');
title(['Vrms Distribution']);
legend(['\sigma=' num2str(std_vrms) ', \mu = ' num2str(mean_vrms)],'location','south');

subplot(2,2,4);
hist(odat(:,2),50);xlabel('Pavg W');ylabel('# Occurences');
title(['Pavg Distribution']);
legend(['\sigma=' num2str(std_pavg) ', \mu = ' num2str(mean_pavg)],'location','south');
add_print_callbacks;
