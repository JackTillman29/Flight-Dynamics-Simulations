function plot( FLT , oversampleFactor, fsf, fsfu)

if(~exist('oversampleFactor','var'))
    oversampleFactor = 10;
end

if(~exist('fsf','var'))
    fsf = 1;
end

if(~exist('fsfu','var'))
    fsfu = 'Hz';
end

FIR_FFT_OVERSAMPLED  = fft(FLT.tapGains,oversampleFactor*(FLT.nTaps+1));
FIR_FFT_OVERSAMPLED_FREQ = FFT_FreqBinCenters(length(FIR_FFT_OVERSAMPLED),FLT.Fs);

figure(1);
subplot(2,1,1);
plot(FLT.tapGains,'.-');
xlabel('Tap #');
ylabel('Tap Gain');
grid on;
axis tight;
title({'Time Domain Convolution FIR Waveform',['Output Delay: ' num2str(FLT.latency) 's']}); 
hold on;
subplot(2,1,2); 
plot(fsf*FIR_FFT_OVERSAMPLED_FREQ,20*log10(abs(FIR_FFT_OVERSAMPLED)));
grid on;
xlabel(['Frequency (' fsfu ')']);
ylabel('Attenuation (20Log_{10}dB)');
title('Frequency Domain Response of FIR Filter'); 
a_ylim = ylim;
hold on;

line(fsf*[FLT.Fs/2 FLT.Fs/2],[a_ylim(1) a_ylim(2)],'Color','k','LineWidth',3);
hold off;


end

