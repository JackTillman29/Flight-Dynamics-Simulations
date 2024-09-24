function x = sawOsc(f,T,dt,phsOffset)
% x = sawOsc(f,T,dt,phsOffset)
%   phsOffset is number between 0 and 1.
if(~exist('phsOffset'))
    phsOffset = 0;
end

Fs = 1/dt;

% a better way is to compute an array of PHASE INCREMENTS (per sample) and 
% then accumulate them using cumsum.
samplesPerCycle = (1./f) / dt;
cyclesPerSample = 1./samplesPerCycle;
phsInc = cyclesPerSample;

N = floor(T/dt);

if(length(phsInc) == 1)
    k = [0:N-1];
    phs = phsInc .* k;
else
    if(length(phsInc) >= N)
        phs = cumsum(phsInc(1:N));
    else
        error(['input f needs to be a scalar or ',num2str(N),' elements long (based on input T and dt)'])
    end
end

x = 2*mod(phs - phsOffset,1) - 1;

end