function x = cosOsc(f,T,dt,phsOffset)
%  x = sinOsc(f,T,dt,phsOffset)

% in this routine, phase offset is ALWAYS POSITIVE, and between 0 and 1
if(~exist('phsOffset'))
    phsOffset = 0;
else
    phsOffset = mod(abs(phsOffset),1);
end

t = [0:dt:(T-dt)];

x = cos(2*pi*(f.*t + phsOffset));

end