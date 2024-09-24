function [bpsk_shifter,bpsk_ampl] = SimplePSK(key,t,varargin)
% help: bpsk_shifter = SimpleLfm(key,t)
% maps to zero to 180 degrees
% use bpsk_signal = complex_carrier_signal*exp(i.*(bpsk_shifter*pi))
% Calculate # of chips based on key length
nChips = length(key);
% Calculate # of chips based on key length
pShift  = round(length(t)/nChips);

bpsk_shifter = 0*t;
bpsk_ampl = 0*t;

iValue = 1;
iRegister = 1;

for kshift = 1:nChips
    iStart = 1+pShift*(kshift-1);
    iStop = iStart + pShift - 1;
    if(iStop > length(t))
        iStop = length(t);
    end
    
    binaryValue = key(kshift);
    bpsk_shifter(iStart:iStop) = binaryValue;
    
    if(nargin > 2)
        bpsk_ampl(iStart:iStop) = varargin{1}(kshift);
    else
        bpsk_ampl(iStart:iStop) = 1;
    end
    
end


end