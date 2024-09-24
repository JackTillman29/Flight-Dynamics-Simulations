close all;
clear all;
clc;

a =load('RawData.mat');
b =load('RawData2.mat');

f_a    = a.RawData(:,1);
attn_a = a.RawData(:,2);

f_b    = b.RawData2(:,1);
attn_b = b.RawData2(:,2);

f_c = 10*log10(logspace(0,2,100));

attn_ac = interp1(f_a,attn_a,f_c);
attn_bc = interp1(f_b,attn_b,f_c);

attn_combined = 10*log10(10.^(0.1*attn_ac) + 10.^(0.1*attn_bc));

figure(1);
loglog(...
    10.^(0.1*f_a),10.^(0.1*attn_a),'sq', ...
    10.^(0.1*f_b),10.^(0.1*attn_b),'sq');
xlabel('Frequency (GHz)');
ylabel('dB/km');
title('Atmospheric Attenuation Table');
hold on;
loglog(...
    10.^(0.1*f_c),10.^(0.1*attn_ac), ...
    10.^(0.1*f_c),10.^(0.1*attn_bc), ...
    10.^(0.1*f_c),10.^(0.1*attn_combined),'LineWidth',2);
legend('H_20','O_2','H_20 (Interp)','O_2 (Interp)','Combined','Location','best');
grid on;
fig_props = PrepForPrint;
PrepForPrint(1,fig_props);
% convert data to GHz,dbPerkm
InputGHz = 10.^(0.1*f_c);
OutputdBPerKm = 10.^(0.1*attn_combined);

inan = find(isnan(OutputdBPerKm));
for k = 1 : length(inan)
    try
        OutputdBPerKm(inan(k)) = OutputdBPerKm(inan(k)+1);
    catch
        OutputdBPerKm(inan(k)) = OutputdBPerKm(inan(k)-1);
    end
end

save('AtmosphericTable.mat','InputGHz','OutputdBPerKm');

% dB per kilometer
lookUpExample = interp1(InputGHz,OutputdBPerKm,10.0)

print -dpng -f1 output_table.png -r180
