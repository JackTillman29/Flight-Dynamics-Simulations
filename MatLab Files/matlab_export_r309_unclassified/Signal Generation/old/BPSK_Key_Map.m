% Code to map out MLS sequences using BPSK keys
close all
addpath('z:\sawmillkd\MATLAB\PulsesAndWaveforms');
nChips     = 1023;
nRegisters = log10(nChips+1) / log10(2);

lengthMap = zeros(nRegisters);

for iReg = 1:nRegisters
    for jReg = 1:nRegisters
        if ( iReg == jReg ) 
            keylength = nan;
        else
            [key,keylength] = BPSK_Key(nChips,[iReg jReg]);
        end
        lengthMap(iReg,jReg) = keylength;
    end
end

disp(lengthMap);

figure(1);
hImage = imagesc(lengthMap);

set(hImage,'AlphaData',~isnan(lengthMap));
set(gca,'YDir','normal');
xlabel('Shift Register Feedback Tap #1');
ylabel('Shift Register Feedback Tap #2');
title('BPSK Shift Register Sequence Length vs Tap Configuration');
colorbar;

for iReg = 1:nRegisters
    for jReg = 1:nRegisters
        text(iReg,jReg,num2str(lengthMap(iReg,jReg),'%d'), ...
            'color','w','HorizontalAlignment','center', ...
            'FontSize',8);
    end
end
