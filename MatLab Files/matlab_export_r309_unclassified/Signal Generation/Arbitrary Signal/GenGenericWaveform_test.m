close all; clear all; clc;

%==========================================================================
% define wideband signal
nFreqs   = 10;
centerRf = 50e6;
spreadRf = 40e6;
duration = 1e-3;

%==========================================================================
% define sub-band of interest
bwOfInterest = 5e6;
rfOfInterest = 50e6;

Fs = 30*bwOfInterest; dt = 1/Fs;
if_filter_freq = linspace(-bwOfInterest/2,bwOfInterest/2,1e3);

% [OPTIONAL] define an IF filter (using butterworth approximation)
butterOrder = 4;
if_filter_gain = ... % magnitude function
    1./(1 + (2*pi*(if_filter_freq./(2*bwOfInterest))).^(2*butterOrder));

% figure;
% plot(if_filter_freq,10*log10(if_filter_gain))
% return

% bwOfInterest = 2e6;
% rfOfInterest = 2.5e9;
% nFreqs   = 10;
% centerRf = 2.5e9;
% spreadRf = 1.5e9;
% duration = 1e-3;

%==========================================================================
% build wideband signal generation tables
freqTable = centerRf + spreadRf/2*(2*rand(1,nFreqs)-1);
timeTable = linspace(0,duration,nFreqs);

timeTableA = timeTable;
ampTable   = ones(1,nFreqs);

%% GENERATE SIGNALS
% ENTIRE SIGNAL
xfull = GenGenericWaveformIF(timeTable,freqTable,timeTableA,ampTable,dt);
Nfull = length(xfull);

% SUB-BAND SIGNAL
x = GenGenericWaveformIF(timeTable,freqTable,timeTableA,ampTable,dt,...
    'rf',rfOfInterest,'bw',bwOfInterest);
N = length(x);

% SUB-BAND SIGNAL (WITH BUTTERWORTH IF FILTER APPLIED)
xbutter = GenGenericWaveformIF(...
    timeTable,freqTable,timeTableA,ampTable,dt,...
    'rf',rfOfInterest,'bw',bwOfInterest,...
    'if_filter_freq',if_filter_freq,'if_filter_gain',if_filter_gain);
Nbutter = length(xbutter);

%% PLOTS
% plot freq vs time (true "waterfall" plot)
figure
plot(freqTable,timeTable,'k')
set(gca,'YDir','reverse')
hold on
ylims = ylim;
plot((rfOfInterest-bwOfInterest/2)*ones(1,2),ylims,'r-')
plot((rfOfInterest+bwOfInterest/2)*ones(1,2),ylims,'r-')

% plot time-domain waveforms
figure
subplot(3,1,1)
plot(real(x))
title('Sub-band signal [w/ Perfect IF Filter]')

subplot(3,1,2)
plot(real(xbutter))
title('Sub-band signal [w/ Butterworth IF Filter)')

subplot(3,1,3)
plot(real(xfull))
title('Full spectrum signal')

% plot each signal's spectrum on single plot
xf       = WaveObj(x,[],Fs,N);
xbutterf = WaveObj(xbutter,[],Fs,Nbutter);
xfullf   = WaveObj(xfull,[],Fs,Nfull);

xf.name       = 'Sub Signal [perfect IF filter]';
xbutterf.name = 'Sub Signal [Butterworth IF filter]';
xfullf.name   = 'Full signal';

plot([xfullf xf xbutterf])

% plot each signal's STFT in single figure
figure('Position',[95 73 1294 678]); add_print_callbacks;
hax_stft(1) = subplot(1,3,1);
hax_stft(2) = subplot(1,3,2);
hax_stft(3) = subplot(1,3,3);

winLen = 1e3;
winFun = [];
overlap = 0.25;
nfft = max([N Nfull Nbutter]);
tSc = 1e6; tU = '\mus';
fSc = 1e-6; fU = 'MHz';
f0 = 0;

clims = [-60 0];

xfullf.stft(winLen,winFun,overlap,nfft,tSc,tU,fSc,fU,f0,hax_stft(1));
title(hax_stft(1),'Full Signal STFT')
set(hax_stft(1),'YDir','reverse'); colorbar off; caxis(clims);

xf.stft(winLen,winFun,overlap,nfft,tSc,tU,fSc,fU,f0,hax_stft(2));
title(hax_stft(2),'Sub Signal [perfect IF filter] STFT')
set(hax_stft(1),'YDir','reverse'); colorbar off; caxis(clims);

xbutterf.stft(winLen,winFun,overlap,nfft,tSc,tU,fSc,fU,f0,hax_stft(3));
title(hax_stft(3),'Sub Signal [Butterworth IF filter] STFT')
set(hax_stft(1),'YDir','reverse'); colorbar off; caxis(clims);

xlims = xlim(hax_stft(3));
ylims = ylim(hax_stft(3));

set(hax_stft(1),'XLim',xlims,'YLim',ylims,'YDir','reverse')
set(hax_stft(2),'XLim',xlims,'YLim',ylims,'YDir','reverse')
set(hax_stft(3),'XLim',xlims,'YLim',ylims,'YDir','reverse')

