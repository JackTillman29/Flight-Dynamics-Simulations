function b = hgain(a)
       b = ((2.0*a^3 - 4.0*a^2 + sqrt(4.0*a^6 - 64.0*a^5 + 64.0*a^4))/(8.0*(1.0 - a)));
end