%a = 0.496;
a = 0.32;
%b = 0.0;


Ts = 0.1;
%b = 6.15e-2;%6.15E-2/Ts;
b = hgain(a);

freqVectorHz = 0:0.01:(0.5/Ts);

% convert frequency to fraction of sample rate * twopi
fracSampleRate = 2*pi*freqVectorHz ./ (1/Ts);

% Build discrete transfer function for alpha-beta
tf.num = [ 0   a+b      -a ];
tf.den = [ 1   a+b-2   1-a ];

tf.zeros = roots(tf.num);
tf.poles = roots(tf.den);

% Build z-domain transfer function (for reference)
zd = exp(1i*fracSampleRate);
output = polyval(tf.num,zd)./polyval(tf.den,zd);

figure(1);
subplot(2,1,1);
semilogx(freqVectorHz,20*log10(abs(output)),'LineWidth',2);
grid on;xlabel('Frequency (Hz)');ylabel('||Y|| /||U|| (dB)');
title({'Open Loop Response',['\alpha=' num2str(a) ', \beta=' num2str(b)]});
%axis tight;
hold on;
subplot(2,1,2);
semilogx(freqVectorHz,180/pi*unwrap(angle(output)),'LineWidth',2);
grid on;xlabel('Frequency (Hz)');ylabel('Phase Lag (deg)');
%axis tight;
hold on;


%% Closed loop
tfcl.num = tf.num;
tfcl.den = [ ...
    tf.num(1)+tf.den(1) ...
    tf.num(2)+tf.den(2) ...
    tf.num(3)+tf.den(3) ];

tfcl.zeros = roots(tfcl.num);
tfcl.poles = roots(tfcl.den);

% Build z-domain transfer function (for reference)
zd = exp(1i*fracSampleRate);
output = polyval(tf.num,zd)./polyval(tf.den,zd);

figure(2);
subplot(2,1,1);
semilogx(freqVectorHz,20*log10(abs(output)),'LineWidth',2);
grid on;xlabel('Frequency (Hz)');ylabel('||Y|| /||U|| (dB)');
title({'Closed Loop Response',['\alpha=' num2str(a) ', \beta=' num2str(b)]});
%axis tight;
hold on;
subplot(2,1,2);
semilogx(freqVectorHz,180/pi*unwrap(angle(output)),'LineWidth',2);
grid on;xlabel('Frequency (Hz)');ylabel('Phase Lag (deg)');
%axis tight;
hold on;


%% Figure

addpath('z:\sawmillkd\MATLAB\Controls');
dstep(tfcl.num,tfcl.den,1000,10);