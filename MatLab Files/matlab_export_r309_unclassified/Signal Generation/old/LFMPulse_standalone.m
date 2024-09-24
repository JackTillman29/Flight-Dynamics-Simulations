function pulse = LFMPulse(f1,f2,pw,ovrsamp)
% pulse = SquarePulse(freqHz,pwSec,ovrsamp=10)

    if ( nargin == 3 )
        ovrsamp = 10.0;
    end
    Fs = ovrsamp * max(f1,f2);
    Ts = 1 / Fs;
    t = 0:Ts:pw;
    %pulse = exp(1i * 2 * pi * f * t);
    %pulse = chirp(t,f1,t(end),f2);
    f = interp1([0 pw],[f1 f2],t);
    pulse = exp(1i * 2 * pi .* f .* t);
end