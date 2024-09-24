Fs = 1000; %Hz
Ts = 1.0 / Fs;

t = 0:Ts:10; %Sec

frequency = 700; %Hz
amplitude = 1;

signal = amplitude * cos(2*pi*frequency*t);

fft_signal = fftshift(fft(signal));
fft_freq = linspace(-Fs/2,Fs/2,length(t));

figure;
plot(fft_freq,(abs(fft_signal))./length(t));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
ylim([-.1 0.6]);
pp=PrepForPrint();
PrepForPrint(get(gcf,'Number'),pp);
set(gcf,'Color','w');

