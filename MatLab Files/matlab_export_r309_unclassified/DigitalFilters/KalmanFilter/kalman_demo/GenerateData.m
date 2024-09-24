function[t,trueData,measData] = GenerateData(xi,yi,xvi,yvi,siteX,siteY,stdRange,stdAngle,dt)
% generates perfect trajectory data for cannon launched projectile
% also measured polar data
tend = 2*yvi/9.80665;

t = [0:dt:tend]';
n = length(t);

trueData.position_x  = xi + xvi*t;
trueData.velocity_x  = xvi * ones(n,1);
trueData.position_y  = yi + yvi*t -9.80665/2 * t.^2;
trueData.velocity_y  = yvi - 9.80665*t;

measData.range = sqrt((trueData.position_x - siteX).^2 + (trueData.position_y - siteY).^2) + ...
    stdRange * randn(n,1);
measData.angle = atan2(trueData.position_y - siteY,trueData.position_x - siteX) + ...
    stdAngle * randn(n,1);
end