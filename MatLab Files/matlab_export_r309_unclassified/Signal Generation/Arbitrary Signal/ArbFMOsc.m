function x = ArbFMOsc(freqArray,dt)
% Arbitrary Frequency Modulation Oscillator
%  freqArray:

% freq = dphs/dt (approx)
%  given freq function (freqArray)
%  --> phs = integrate(f*dt)
phs = 2*pi*cumsum(freqArray*dt);  % phase accumulator
x = exp(1i*phs);

end