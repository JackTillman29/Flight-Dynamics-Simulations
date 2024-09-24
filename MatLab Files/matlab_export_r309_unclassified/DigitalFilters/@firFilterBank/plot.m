function plot( FLTB , oversampleFactor )

color_plot = [ ...
    1 0 0;
    0 1 0;
    0 0 1 ];

if(~exist('oversampleFactor','var'))
    oversampleFactor = 10;
end
Fs = get(FLTB.filters{1},'Fs');
figure;
for k = 1 : length(FLTB.filters)
    FIR_FFT_OVERSAMPLED         = fft( ...
        get(FLTB.filters{k},'tapGains'), ...
        oversampleFactor*(get(FLTB.filters{k},'nTaps')+1) ...
        );
    FIR_FFT_OVERSAMPLED_FREQ    = FFT_FreqBinCenters( ...
        length(FIR_FFT_OVERSAMPLED), ...
        Fs ...
        );
    hold on;
    plot(FIR_FFT_OVERSAMPLED_FREQ,10*log10(abs(FIR_FFT_OVERSAMPLED)),'Color', ...
        color_plot(mod(k,3)+1,:));
end




grid on;
xlabel('Frequency (Hz)');
ylabel('Attenuation (10Log_{10}dB)');
title('Frequency Domain Response of FIR Filter Bank'); 
a_ylim = ylim;
hold on;

line([Fs/2 Fs/2],[a_ylim(1) a_ylim(2)],'Color','k','LineWidth',3);
hold off;

end

