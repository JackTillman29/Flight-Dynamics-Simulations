% initial conditions
m = 10;
vx = 100;
vy = 100;
t = 0:0.01:20;
g = -9.80665;
radarx = 1;
radary = 0;

% matrices
F = [ ...
    0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 0];

B = [...
    0;
    0;
    0;
    g];

C = eye(4);

D = 0;

sys1 = ss(F,B,C,D);

dsys1 = c2d(sys1,0.01,'tustin');


y = zeros(length(t),4);
x = [0 vx 0 vy]';
y(1,:) = x';
x = inv(dsys1.c) * x;
for k = 2:length(t)
    x = dsys1.a * x + dsys1.b;
    y(k,:) = dsys1.c*x+dsys1.d;
end

xm = y(:,1) - radarx;
ym = y(:,3);

rm = sqrt((xm.^2+ym.^2)) + 50*randn(length(t),1);
am = atan2(ym,xm) + 1.0*pi/180*randn(length(t),1);

xm_noise = rm.*cos(am);
ym_noise = rm.*sin(am);