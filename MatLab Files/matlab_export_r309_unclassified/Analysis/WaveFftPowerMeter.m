function [pm,pmf] = WaveFftPowerMeter(cpsd,frq,bufw)
% function [pm,pmf] = WaveFftPowerMeter(cpsd,frq,bufw)
% returns the "sectorized" cumulative power spectral density
% this better informs the power in a given band.
    pmf = frq(bufw:end)  - (bufw/2-1)*(frq(2)-frq(1));
    pm  = cpsd(bufw:end) - cpsd(1:(end-bufw+1));
end