function y = LogAmp(x, xIntercept_dBm, outSlope_mV_per_dB, dynamicRange_dB, varargin)
% implements a base-10 logarithmic amplifier
%   y = LogAmp(x, xIntercept_dBm, outSlope_mV_per_dB, dynamicRange_dB, varargin)
% INPUTS:
%   x:                  input signal (can be a scalar, vector, or matrix)
%   xIntercept_dBm:     the input signal strength that produces an output
%                       of zero volts
%   outSlope_mV_per_dB: slope of the 
%   dynamicRange_dB:    
% define intercept voltage
Vx_dBm       = xIntercept_dBm;     % dBm input
Vy_mV_per_dB = outSlope_mV_per_dB; % mV/dB
dynRng_dB    = dynamicRange_dB;    % dB

Vy_mV_per_decade = Vy_mV_per_dB.^2;

Vx1 = 10^((Vx_dBm-30)/20);
Vx2 = 10^((Vx_dBm+dynRng_dB-30)/20);

y = Vy_mV_per_decade*1e-3 .* log10(x ./ Vx1);

end