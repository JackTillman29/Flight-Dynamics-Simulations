function outputStructure = WaveSTFTStruct(tdSignal,tdWindow,Fs,shortTime,shortOverlapPct)
%                              WaveSTFTStruct(tdSignal,tdWindow,Fs,shortTime,shortOverlapPct)
% function [outputStructure] = WaveSTFTStruct(ys,@hamming,1000,.25,.25);
% Short Time Fourier Transform (generates spectragraph)
% no output arguments make a plot.




Ts = 1/Fs;
t = Ts * ((1:length(tdSignal)) - 1);


tWindow = shortTime;
pctOvlp = shortOverlapPct;
nWindow = tWindow ./ Ts-1;
nOffset = ceil(nWindow*(1-pctOvlp));

% poor mans way to get how many overlapping windows
for idx = 1:length(t)
    if(idx == 1)
        start = 1;
    else
        start = stop + 1 - nOffset;
    end
    stop  = start + nWindow;
    if(stop > length(t))
        break;
    end
end
nWindows = idx - 1;
% form image data
outputStructure.fftdata = complex(zeros(nWindows,nWindow+1));
outputStructure.psddata = complex(zeros(nWindows,nWindow+1));
outputStructure.timedata = zeros(nWindows,1);

% now compute the ffts
for iwin = 1:nWindows
    if(iwin == 1)
        start = 1;
    else
        start = stop + 1 - nOffset;
    end
    stop  = start + nWindow;
    midp = floor(start + (stop-start)/2);
    outputStructure.timedata(iwin) = t(midp);
    
    
    [sig_fft, sig_fft_freq, sig_psd] = WaveFft(tdSignal(start:stop),tdWindow(length(start:stop))',Fs);
    outputStructure.fftdata(iwin,:) = sig_fft;
    outputStructure.psddata(iwin,:) = sig_psd;
    if(iwin==1)
        outputStructure.freqdata = sig_fft_freq;
    end
    
end

if(nargout == 0)
    figure;
    imagesc(outputStructure.freqdata,outputStructure.timedata,10*log10(outputStructure.psddata));
    colorbar;
    xlabel('Frequency (Hz)');
    ylabel('Time (s)');
    title({['Short Time PSD'],['Window Size: ' num2str(shortTime) 's, % Overlap: ' num2str(100*shortOverlapPct) ' %']});
end

end