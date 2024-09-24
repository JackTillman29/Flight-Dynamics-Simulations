close all; clear all; clc

r2d = 180/pi;
d2r = 1/r2d;

rng = 311536; % meters
ang = -2.0;    % deg
bw  = 4.2;    % deg

rng = 311536
ang = 2.1;
bw  = 4.2;

alt = 9144;
beamwidth_rad = 1 * d2r;

N = 100;
gnd = linspace(50e3,10e3,N); % downrange

rngToGnd = sqrt(alt.^2 + gnd.^2);
angToGnd = atan2(alt,gnd);

figure;

[vbwGndRng,hbwGndRng,data] = ...
    ProjectionOnGround(rngToGnd,angToGnd,beamwidth_rad);

figure;
subplot(3,1,1)
plot(rngToGnd * 1e-3); title('Range to Ground along LOS [km]')
subplot(3,1,2)
plot(angToGnd * r2d); title('Angle to Ground along LOS [deg]')
subplot(3,1,3)
% plot(bwGndRng); title('Range Width of Ground-Projected BW')
plot([data.vbwRng1.' data.vbwRng2.']); title('Range Width of Ground-Projected BW'); legend('Inner V-proj','Outer V-proj')


figure;
plot(data.hbwRng)
hold on;
plot(data.vbwRng)
legend('H-proj','V-proj')