close all;
clear all;
clc;
addpath('I:\MATLAB\DSP');
Fif = 100e6;
Fs  = 8*Fif;

dt = 1/Fs;

pw = 128e-6;
chirp = 40e6;

[TxRef,TimeRef,FreqRef] = SimpleLfm(Fif-0.5*chirp,Fif+0.5*chirp,pw,dt);
TxRef = single(TxRef);
TimeRef = single(TimeRef);
FreqRef = single(FreqRef);

plat2satDoppler = 6e3;
plat2gndDoppler = 0*9.8e3;

RxDirect = TxRef .* exp(1i*2*pi*plat2satDoppler*TimeRef);




RxTwoWay = TxRef .* exp(1i*2*pi*plat2gndDoppler*TimeRef);

filt = PulseCompressionXCORRsStruct(RxTwoWay,RxDirect,'no_align');

figure;
plot(1e6*dt*filt.kleadlag,20*log10(abs(filt.signalCrossCor)));
T2 = (filt.kleadlag - filt.kleadlag(1)) .* dt;
RxDirectLO =  exp(1i*2*pi*(plat2satDoppler+Fif)*T2);
downc = filt.signalCrossCor ./ RxDirectLO;
%%
[sig_fft, sig_fft_freq, sig_psd, sig_cpsd] = WaveFft(downc,0*downc+1,Fs,'onesided');
figure;
subplot(2,1,1);
plot(1e6*dt*filt.kleadlag,real(downc));
subplot(2,1,2);
plot(sig_fft_freq,abs(sig_fft));
