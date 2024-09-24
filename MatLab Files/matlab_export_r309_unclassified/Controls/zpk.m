function [z,p,k] = zpk(num,den)
% This function decomposes a polynomial transfer function into
% pole/zero/gain form.
% 
% k*(s-z1)(s-z2)...(s-zn)
% -------------------------------
%   (s-p1)(s-p2)...(s-pn)

% K. Sawmiller, 2011
%
% Syntax: [z,p,k] = zpk(num,den)

% Compute the roots of the numerator and denominator
    z = roots(num);
    p = roots(den);

% Reform the polynomial based upon the roots computed above
% This is needed to figure out the "lost" polynomial gain
    num_poly = poly(z);
    den_poly = poly(p);
    
% Compute the ratio of the new polynomial to the original in order
% to get the gain factors for both numerator and denominator
    [dummy1,dummy2,k_num] = residue(num,num_poly);
    [dummy1,dummy2,k_den] = residue(den,den_poly);

% Combine the gains into a single value applied to the numerator
    k = k_num / k_den;
    
end