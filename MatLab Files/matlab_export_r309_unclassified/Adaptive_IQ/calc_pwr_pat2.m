function [pwr,angdeg] = calc_pwr_pat2(s,D,lambda)
    ang = linspace(-pi/2,pi/2,8191);
    N = length(s);
    nVec = 0:N-1;
    y = zeros(1,length(ang));
    for k = 1 : length(y)
        y(k) = s'*exp(1i*2*pi*nVec*D*sin(ang(k))./lambda).';
    end
    pwr = 20*log10(abs(y));
    angdeg = 180/pi*ang;
end