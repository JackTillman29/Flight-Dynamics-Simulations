% Plot to the nyquist
omega = 0.001:0.005:pi;
N = 1;
Ts = 4e-3;
Fs = 1.0 / Ts;
pcolor = 'k';
imResult = (1/N)*(1-exp(-1i*omega*N))./(1-exp(-1i*omega));
mag = abs(imResult);
phase = unwrap(angle(imResult));

% -- Plot Magnitude -- %
subplot(2,1,1);
hold on;
semilogx(1e-3*omega*Fs/(2*pi),20*log10(mag),pcolor);
grid on;axis tight;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title({['Bode Plot of ' num2str(N) ' Sample Averaging Filter']});
hold off;

% -- Plot Phase -- %
subplot(2,1,2);
hold on;
semilogx(1e-3*omega*Fs/(2*pi),unwrap(phase*180/pi),pcolor);
grid on;axis tight;
xlabel({'Frequency (kHz)'});
ylabel('Phase (Deg)');
hold off;
