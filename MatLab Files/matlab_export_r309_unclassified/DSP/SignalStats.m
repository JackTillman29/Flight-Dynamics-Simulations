function signalStats = SignalStats(x,Fs)

%% compute time energy density function
%   defined: |s(t)|.^2
xtd = abs(x).^2; % time density function
E = sum(xtd)*dt
xtd = xtd / E; % divide by total signal energy to normalize density function

%% time stats (mean, var, std)
meanTint = cumsum(t.*xtd.*dt);     % <w>
ip_tsq   = meanTint.^2;            % <w>.^2
ip_t2    = cumsum(t.^2.*xtd.*dt);  % <w^2>  ip == inner product

% variance time == <w^2> - <w>^2
varTint = ip_t2 - ip_tsq;
% std time (RMS Time Duration)
stdTint = sqrt(varTint);

signalStats.timeAvg = meanTint(end);
signalStats.timeVar = varTint(end);
signalStats.timeStd = stdTint(end);

%% compute frequency energy density function
%   defined: |S(f)|.^2
ws = fftshift(fft(s,N)/N);           % fourier transform of signal
df = Fs / N; f = [-Fs/2:df:Fs/2-df]; % frequency vector [-Fs/2,Fs/2)
dw = df*2*pi/Fs; w = f*2*pi/Fs;      % rad/sec vector   [-pi  ,pi)

signalStats.eds  = abs(ws).^2;     % energy density spectrum
signalStats.psd  = eds / df;       % power spectral density
signalStats.ceds = cumsum(eds);    % cumulative energy spectra
signalStats.cpsd = cumsum(psd)*df; % cumulative power spectra
signalStats.sigPower = cpsd(end)   % signal power computed from integrating under all freq in PSD

%% freq stats (mean, var, std)
meanWint = cumsum(w.*signalStats.eds.*dt);     % <w>
ip_wsq   = meanWint.^2;                        % <w>.^2
ip_w2    = cumsum(t.^2.*signalStats.eds.*dt);  % <w^2>  ip == inner product

% variance time == <w^2> - <w>^2
varWint = ip_w2 - ip_wsq;
% std time (RMS Time Duration)
stdWint = sqrt(varTint);

signalStats.freqAvg = meanWint(end) / (2*pi);
signalStats.freqVar = varWint(end) / (2*pi);
signalStats.freqStd = stdWint(end) / (2*pi);
signalStats.BW_RMS = signalStats.freqStd;

%% RMS Bandwidth


end