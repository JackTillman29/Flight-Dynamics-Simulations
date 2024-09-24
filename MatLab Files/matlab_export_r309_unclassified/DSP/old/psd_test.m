close all;clear all;clc;
dt = 0.0025;
t  = 0:dt:(50-dt);



F1 = 100;
F2 = 300;

% Total Power
PT = 1;

% True PSD, W/Hz
TRUE_PSD = PT; % / (F2-F1);


% True ASD, V/Hz
%TRUE_ASD = sqrt(TRUE_PSD * DF)


NFFT = length(t);
Fs   = 1.0 / dt;
DF   = Fs / NFFT
PSDSF= 1.0 / DF
FV=Fs/NFFT*(0:(NFFT-1));

% True A/Bin
%TRUE_APB = TRUE_ASD * DF
TRUE_PSD2 = TRUE_PSD * DF / (F2-F1);


fy = 0*(0:(NFFT - 1));
fBounds = find(FV >= F1 & FV <= F2);
fy(fBounds) = sqrt(TRUE_PSD2) * NFFT;

y = ifft(fy);

fy = fy./NFFT;


subplot(3,1,1);
plot(Fs/NFFT*(0:(NFFT-1)),abs(fy),'.-');
subplot(3,1,2);
plot(Fs/NFFT*(0:(NFFT-1)),fy.*conj(fy)./DF,'.-');
subplot(3,1,3);
plot(Fs/NFFT*(0:(NFFT-1)),cumsum(fy.*conj(fy)),'.-');