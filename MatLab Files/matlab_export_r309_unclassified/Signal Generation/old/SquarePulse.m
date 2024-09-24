function pulse = SquarePulse(f,pw,ovrsamp)
% pulse = SquarePulse(freqHz,pwSec,ovrsamp=10)

    if ( nargin == 2 )
        ovrsamp = 10.0;
    end
    Fs = ovrsamp * f;
    Ts = 1 / Fs;
    t = 0:Ts:pw;
    pulse = exp(1i * 2 * pi * f * t);
end