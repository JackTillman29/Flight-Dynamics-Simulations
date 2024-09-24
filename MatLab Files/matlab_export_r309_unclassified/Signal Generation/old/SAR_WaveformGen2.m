close all;clear all;clc;
addpath('I:\MATLAB\Windows');

% build baseline transmitted waveform
Fc = 100e6;  % 
Fm = 8 * Fc; % Model Frequency, kHz

dt = (1/Fm);
pw = 128e-6;
chirp = 20e6;
[pulseOut,t] = SimpleLfm(Fc-0.5*chirp,Fc+0.5*chirp,pw,dt);

%% Code to generate full range doppler matrix for a single signal
Fd = 1e3;
FdRef = 1.0e3;
pri_times = 220e-6*(0:31); % start times of each pulse of the 32 pulse sequence
phaseOffsets = single(exp(1i*2*pi*Fd*pri_times));
phaseOffsetsRef = single(exp(1i*2*pi*FdRef*pri_times));

AD_SAMPLES = single(complex(zeros(16,32)));
FFT_SAMPLES = single(complex(zeros(16,32)));

% Determine the timing of the range gate for all 16 pulses

tGateSpacing = round(2*100/3e8/dt);
[inputCrossCorr,kleadlag,kernelAutoCorr,fftInput,fftKernel] = PulseCompressionXCORRs(pulseOut,pulseOut,dt);
mdl = (length(inputCrossCorr)-1)/2+1;
tGate = (-7:8).*tGateSpacing + mdl;
pickoffs = complex(single(zeros(length(phaseOffsets),16)));

for k = 1:length(phaseOffsets)
    DopMod = Fd;
    % Modulate the signal to simulate doppler on the return
    pulseReturn = phaseOffsets(k) .* pulseOut;
    
    % Apply pulse compression
    [inputCrossCorr,kleadlag,kernelAutoCorr,fftInput,fftKernel] = PulseCompressionXCORRs(pulseReturn,pulseOut.*phaseOffsetsRef(k),dt);
    if(~exist('COHO_IF'))
        ttemp = (kleadlag-kleadlag(1))*dt;
        COHO_IF = single(exp(1i*2*pi*Fc*ttemp));
    end
    
    % IQ Demodulation
    iqDemod = (inputCrossCorr' ./ 0.3)./COHO_IF';
    
    for r = 1:16
        leftEdge = tGate(r)-floor(tGateSpacing/2);
        rightEdge = tGate(r)+floor(tGateSpacing/2);
        adReal = mean(real(iqDemod(leftEdge:rightEdge)));
        adImag = mean(imag(iqDemod(leftEdge:rightEdge)));
        AD_SAMPLES(r,k) = single(complex(adReal,adImag));
    end
end

for r = 1:16
    FFT_SAMPLES(r,:) = fft(CosSquared(32).*AD_SAMPLES(r,:))./32;
end

figure;
hs=surf(double(20*log10(abs(FFT_SAMPLES))));
set(hs,'LineStyle',':');
view(0,90);
axis tight;
title('\fontsize{12}\bf{}SA-15 TER Range Doppler Matrix (dBW)');
    ylabel('\fontsize{12}\bf{}Range ADC Gate');
    xlabel('\fontsize{12}\bf{}Doppler FFT Bin');
colorbar;
%%
hr = [];
figure;
set(gcf,'DoubleBuffer','on','BackingStore','on','HitTest','off','Clipping','off','Interruptible','off');
dh = 1/16;
for q = 16:-1:1
    subplot(16,1,q);
    hr = [hr gca];
    plot([real(CosSquared(32).*AD_SAMPLES(q,:))' imag(CosSquared(32).*AD_SAMPLES(q,:))'],'sq-','LineWidth',2,'MarkerSize',4);
end

for q = 1:16
    set(hr(q),'YTick',[],'XTick',[], ...
        'Position',[0 (q-1)*dh 1 dh], ...
        'Box','on','HitTest','off','Interruptible','off','Clipping','off');
    axes(hr(q));
    axis tight;
    ylim([-1 1]);
end

