close all;
clear;
clc;

addpath('Z:\File Cabinet\MATLAB\DSP');

% Range in which the waveform will cover
RCPI = 65000;

% Compute waveform parameters
CPIDUR = 1e-6 * RCPI / 150;
FS = 2.4e6;
LFMBW = 2.4e5;
LFMPW = 64e-6;
CPI_NSAMPLES = floor(CPIDUR * FS);
NoiseFig = 4.5;%dB
NoisePow = 1*1.38e-23 * LFMBW * 290 * 10^(NoiseFig/10);


% Create I/Q reference waveform to be transmitted
[txWaveform]=SimpleLfm( ...
    -LFMBW/2, ...
    LFMBW/2, ...
    LFMPW, ...
    1/FS);

PULSE_NSAMPLES = length(txWaveform);


% Create the return CPI as thermal noise for now
rxSignal = sqrt(NoisePow) * exp(1i*2*pi*rand(1,CPI_NSAMPLES));

% Drop some target returns in the CPI
tgt1 = 400;
tgt2 = 450;
tgt1_SNR = 10^(25/10);
tgt2_SNR = 10^(10/10);
rxSignal(tgt1:(tgt1-1+PULSE_NSAMPLES)) = ...
    rxSignal(tgt1:(tgt1-1+PULSE_NSAMPLES)) + ...
    sqrt(NoisePow*tgt1_SNR)*txWaveform;
rxSignal(tgt2:(tgt2-1+PULSE_NSAMPLES)) = ...
    rxSignal(tgt2:(tgt2-1+PULSE_NSAMPLES)) + ...
    sqrt(NoisePow*tgt2_SNR)*txWaveform;

% Perform basic matched filtering (for reference)
matched = PulseCompressionXCORRsStruct( ...
    rxSignal, txWaveform, 'no_align');


% Create the convolution matrix
C = convmatrix(txWaveform,rxSignal,'trimmed');

% Perform the Singular Value Decomposition (SVD)
[U,S,V] = svd(C);

record_video = 1;

%%

nSigmas = 1:length(diag(S));
first_run = 1;
for nSigma = nSigmas
    
    %InvC = pinv_svd(C,nSigma);
    InvC = pinv_svd_ext(U,S,V,nSigma);
    LSPC = abs(InvC * rxSignal');
    
    if(first_run)
        figure;
        subplot(1,2,1);
        plot(10*log10(diag(S)));
        hL = line([nSigma nSigma],ylim,'Color','r');
        ylim([-60 20]);
        grid on;
        title('Convolution Matrix Singular Values');
        ylabel('10log_{10} dB');
        
        subplot(1,2,2);
        
        hP = plot(10*log10([abs(LSPC) abs(matched.signalCrossCor(77+3:end-77+3).')]));
        hT = title(['Using ' num2str(nSigma) ' of ' num2str(length(C)) ' \sigma']);
        grid on;
        set(hP(2),'Color',0.75*[1 1 1]);
        xlim([200 700]);
        ylim([-95 -60]);
        drawnow;
        first_run = 0;
        pause;
        if(record_video)
            vidobj = VideoWriter('output.avi','MPEG-4');
            open(vidobj);
            vidobj.writeVideo(getframe(gcf));
        end
    else
        set(hP(1),'YData',10*log10(abs(LSPC)));
        set(hT,'String',['Using ' num2str(nSigma) ' of ' num2str(length(C)) ' \sigma']);
        set(hL,'XData',[nSigma nSigma]);
        drawnow;
        if(record_video)
            vidobj.writeVideo(getframe(gcf));
        end
    end
end
if(record_video)
close(vidobj);
    end
