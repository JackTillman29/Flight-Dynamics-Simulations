function dBPower = GetPowerAtRange(rcs,pt,gt,lt,R,freq)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dBPower = pt + gt + lt - 20*log10(R) - 10*log10(4*pi) + rcs - 20*log10(R) - 10*log10(4*pi) + 20*log10(3e8/freq) - 10*log10(4*pi) + gt;



end

