function FLTB = firFilterBank(Fs,nChannels,nTaps,windowString,halfOption)
% function FLTB = firFilterBank(Fs,nChannels,nTaps,windowString,halfOption)
%                 Fs:           Digital Sample Rate
%                 nChannels:    # of channels (full sample rate)
%                 nTaps:        Tap size of FIR filter
%                 windowString: Window applied to FIR filter
%                 halfOption:   1 = zero to Nyquist, 0 = zero to Fs
if(~exist('halfOption','var'))
    halfOption = 0;
end

    
    
    % bank width
    Fbw = Fs / nChannels;
    
    if(halfOption==1)
        nChannels = floor(nChannels/2);
    end
    
    FLTB.filters = cell(nChannels,1);
    
    % build banks
    for k = 1:nChannels
        Fc = Fbw/2.0 + (k-1)*Fbw; % center
        FcStart = Fc - Fbw/2.0;
        FcStop  = Fc + Fbw/2.0;
        %firFilter(fStartHz,fStopHz,Fs,nTaps,windowString)
        disp([FcStart FcStop]);
        FLTB.filters{k} = firFilter(FcStart,FcStop,Fs,nTaps,windowString);
    end
    FLTB.nChannels = nChannels;
    FLTB = class(FLTB,'firFilterBank');
end