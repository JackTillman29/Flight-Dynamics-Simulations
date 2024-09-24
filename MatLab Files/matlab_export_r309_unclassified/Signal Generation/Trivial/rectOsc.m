function x = rectOsc(f,T,dutyCycle,dt,phsOffset)
% will create a DC-compensated pulse waveform
%   "DC-compensated" means that the mean value of one complete cycle is zero.
% ***WARNING*** NOT USEFUL FOR A GATE SIGNAL ***WARNING***
%               --> use clockOsc instead

if(~exist('phsOffset'))
    phsOffset = 0;
end

x1 = sawOsc(f,T,dt,dutyCycle+phsOffset);
x2 = sawOsc(f,T,dt,phsOffset);
x = x1 - x2;

end