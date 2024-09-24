function [b,a] = gen_allpass(norm_freq_ang,norm_freq_mag)
% generates an allpass filter placing poles at specified angles and
% magnitude.
% Syntax:
%   [b,a] = gen_allpass(norm_freq_ang,norm_freq_mag)
%         norm_freq_ang: number between 0 and 1 representing the frequency the
%                        pole-zero pair is positioned for the allpass filter
%         norm_freq_mag:  the strength of this pole-zero pair between 0 and 1
%
% author: Jeff Hole (Booz Allen Hamilton) 2019-3

% norm_freq_ang = rand(1,N);
% norm_freq_mag = 0.8*rand(1,N);

ps = [norm_freq_mag .* exp(-1i*pi*norm_freq_ang)];

% conjugate polynomial and append (ensures REAL roots)
ps = [ps conj(ps)];
ps_poly = [1 -ps(1)];
for k = 2:length(ps)
    ps_poly = conv(ps_poly,[1 -ps(k)]);
end

ps_poly = real(ps_poly);

a = ps_poly;
b = a(end:-1:1);

% numer_roots = roots(b);
% denom_roots = roots(a);




end