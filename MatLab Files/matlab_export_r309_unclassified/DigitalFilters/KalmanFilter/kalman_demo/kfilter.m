% kalman filter

% H matrix
H = @(x,y,xr,yr) [ (x-xr)/sqrt(x^2+y^2) 0 (y-yr)/sqrt(x^2+y^2) 0; ...
    -(y-yr)/(x^2+y^2) 0 (x-xr)/(x^2+y^2) 0];

xhat = [-radarx 100 0 100]';

Pk = [1e-12 0 0 0;
      0 1e2 0 0;
      0 0 1e-12 0;
      0 0 0 1e2];
  
Qk = [
    0 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];
Phi = expm(F*0.01);
I = eye(4);
Rk = [ ...
    50^2 0;
    0 (1*pi/180)^2 ];


kcov = zeros(length(rm),4);
kest = kcov;

kest(1,:) = xhat';
for k = 2:length(rm)
        HH = H(xhat(1),xhat(3),radarx,radary);
        % Solve Riccati Equation for Mk (pre-measurement covariance)
        Mk = Phi * Pk * Phi' + Qk;

        % Solve Riccati Equation for Kk (Kalman Gains)
        Kk = Mk * HH' * inv(HH * Mk * HH' + Rk);

        % Solve Riccati Equation for Pk (post-measurement covariance)
        Pk = (I - Kk * HH) * Mk;
        
        xhat = xhat + Kk * ([rm(k) am(k)]' - [norm(xhat([1 3])) atan2(xhat(3),xhat(1))]') + ...
            g*[0 0 0.01^2/2 0.01]';
        kest(k,:) = xhat';
end