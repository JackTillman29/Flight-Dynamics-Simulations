function f = firFilter(fStartHz,fStopHz,Fs,nTaps,windowString,bstopflag)
%addpath('../Windows');
%addpath('../DSP');
% K. Sawmiller

% Defaults
if(~exist('fStartHz','var'))
    FLT.fStartHz = 200.0;
else
    FLT.fStartHz = fStartHz;
end
if(~exist('fStopHz','var'))
    FLT.fStopHz  = 300.0;
else
    FLT.fStopHz = fStopHz;
end
if(~exist('Fs','var'))
    FLT.Fs  = 1000.0;
else
    FLT.Fs = Fs;
end
if(~exist('nTaps','var'))
    FLT.nTaps = 127;
else
    FLT.nTaps = nTaps;
end
if(~exist('windowString','var'))
    FLT.windowString = 'rectwin(x)';
else
    FLT.windowString = windowString;
end

if(~exist('bstopflag','var'))
    FLT.bstopflag = 0;
else
    FLT.bstopflag = bstopflag;
end

FLT.windowString = strrep(FLT.windowString,'x',num2str(FLT.nTaps));




% Compute constants
nBins = FLT.nTaps + 1;
dt    = 1.0 / FLT.Fs;

% determine the frequency bin centers for this many taps
% FFT_FRQ    = FFT_FreqBinCenters(nBins,FLT.Fs);
df = FLT.Fs / nBins;
FFT_FRQ = [0:(nBins-1)]*df - FLT.Fs/2;

% latency to fill the taps for first output
FLT.latency = FLT.nTaps * dt;


FIR_FFT = zeros(1,nBins);

% this method seems to not give symmetric responses
%[minval,idxstart]=min(abs(FFT_FRQ - FLT.fStartHz));
%[minval,idxstop ]=min(abs(FFT_FRQ - FLT.fStopHz ));

idxstart = find(FFT_FRQ>=FLT.fStartHz,1,'first');
idxstop  = find(FFT_FRQ<=FLT.fStopHz,1,'last');
if(idxstop < idxstart)
    idxstop = idxstart;
end


%disp([idxstart idxstop]);


if(idxstop > nBins/2)
    %error('Pass band extends beyond Nyquist! Check FFT & FIR Relationships!');
    disp('Pass band extends beyond Nyquist! Check FFT & FIR Relationships!');
    idxstop = floor(nBins/2);
end

% set the FIR FFT filter to one in the passband
FIR_FFT(idxstart:idxstop) = 1;
%nBins_ovr_2 = floor(nBins/2);

if(FLT.bstopflag)
    FIR_FFT = real(~FIR_FFT);
end


FIR_FFT(end:-1:(floor(nBins / 2)+2)) = FIR_FFT(2:(ceil(nBins / 2)));
FIR_IFFT = ifft(FIR_FFT);
% throw away point for symmetry
idxthrow = round(nBins / 2 + 1);

FLT.windowFcn = eval(FLT.windowString);

FLT.tapGains = FLT.windowFcn .* real(FIR_IFFT([(idxthrow+1):end 1:(idxthrow-1)]));

f = class(FLT,'firFilter');
end