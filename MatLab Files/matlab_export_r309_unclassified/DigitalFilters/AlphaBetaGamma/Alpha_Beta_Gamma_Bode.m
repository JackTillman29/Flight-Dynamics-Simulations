%%
clear all;
close all;

T=0.025;
alpha=0.4;
beta=0.05;
gamma=0.0;
%dt = .1;
%alpha = 0.32;
%beta = 6.15e-2/dt;
%gamma = 2*2.9494e-3/(dt*dt);

count=0;
for w = 0.0:0.1:pi/T
    count = count + 1;
    
    freq(count) = w / (2 * pi);
    
    Z = exp(i*w*T);
    delta = (Z-1)^3 + alpha*(Z-1)^2 + Z*(Z-1)*(beta + gamma*0.5) + gamma*Z;
    num = Z*(alpha * (Z-1)^2 + (beta + gamma*0.5) * (Z-1) + gamma);
    answer = num / delta;
    mag(count) = abs(answer);
    phase(count) = angle(answer);
end

%SetPlotDefaults(12);

%edata = load('bode_rng.txt');

semilogx(freq, 20*log10(mag),'LineWidth',3);
grid on;
%hold on;
%semilogx(edata(:,1),20*log10(edata(:,3)),'+r');
xlabel('Input Frequency (Hz)');
ylabel('Magnitude (dB, 20 Log10)');
title('Bode Plot, Magnitude');

figure(2);
semilogx(freq, phase*180/pi, 'LineWidth',3);
xlabel('Input Frequency (Hz)');
ylabel('Output/Input Phase (deg)');
title('Bode Plot, Phase');

%figure(3);
%polar(phase,mag);
%grid on;
