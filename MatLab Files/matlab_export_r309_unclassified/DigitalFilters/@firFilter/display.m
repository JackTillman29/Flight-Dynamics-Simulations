function display( FLT )
    disp('FIR Filter Object');
    disp(sprintf('Sample Rate (Hz)      : %f',FLT.Fs));
    disp(sprintf('Start of Passband (Hz): %f',FLT.fStartHz));
    disp(sprintf('End of Passband (Hz)  : %f',FLT.fStopHz));
    disp(sprintf('# of Taps             : %d',FLT.nTaps));
    disp(sprintf('Window Function       : %s',FLT.windowString));
    disp(sprintf('Latency (s)           : %f',FLT.latency));
end