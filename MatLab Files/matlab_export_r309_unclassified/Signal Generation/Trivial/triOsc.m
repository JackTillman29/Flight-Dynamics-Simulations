function x = triOsc(f,T,dt,phsOffset)
%  x = triangle(f,a,dt)

% in this routine, phase offset is ALWAYS POSITIVE, and between 0 and 1
if(~exist('phsOffset'))
    phsOffset = 0;
else
    phsOffset = mod(abs(phsOffset),1);
end

% integrating a DC-compensated rectangular wave gives a triangle wave
dutyCycle = 0.5; % TODO: figure out how to create a DC-compensated triangle
                 %       oscillator for any duty cycle input.
xsquare = rectOsc(f,T,dutyCycle,dt,phsOffset);

samplesPerPeriod = (1./f) / dt;
nRiseSamples = samplesPerPeriod / 2;

% if there is a phase offset, then an initial state must be calculated to
% ensure the triangle waveform is DC-compensated
if(phsOffset < 0.5)
    % neg swing of square wave
    x0 = 4*phsOffset;
else
    % pos swing of square wave
    x0 = 4*(1-phsOffset);
end

x = 2 * cumsum(xsquare./nRiseSamples) - 1 + x0;

end