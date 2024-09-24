function[t,trueData,measData] = GenerateData2(xi,yi,xvi,yvi,siteX,siteY,stdRange,stdAngle,dt,dragc)
% generates perfect trajectory data for cannon launched projectile
% also measured polar data
tend = 2*yvi/9.80665;

t = [0:dt:tend]';
n = length(t);

% just for sizing
    trueData.position_x  = xi + xvi*t;
    trueData.velocity_x  = xvi * ones(n,1);
    trueData.position_y  = yi + yvi*t -9.80665/2 * t.^2;
    trueData.velocity_y  = yvi - 9.80665*t;

    trueData.position_x(1) = xi;
    trueData.position_y(1) = yi;
    trueData.velocity_x(1) = xvi;
    trueData.velocity_y(1) = yvi;
    
for k = 2:length(t)
   % fpa  = atan2(trueData.velocity_y(k-1),trueData.velocity_x(k-1));  % flight path angle
    drag = dragc * sqrt(trueData.velocity_x(k-1)^2+trueData.velocity_y(k-1)^2); % acceleration from drag;
    dragx = -trueData.velocity_x(k-1) * drag;
    dragy = -trueData.velocity_y(k-1) * drag - 9.80665;
    trueData.velocity_x(k) = trueData.velocity_x(k-1) + dragx * dt;
    trueData.velocity_y(k) = trueData.velocity_y(k-1) + dragy * dt;
    trueData.position_x(k) = trueData.position_x(k-1) + trueData.velocity_x(k-1)*dt + dragx * dt^2/2;
    trueData.position_y(k) = trueData.position_y(k-1) + trueData.velocity_y(k-1)*dt + dragy * dt^2/2;
end

if(dragc ~= 0)  % might need to recheck length of data
    row = max(find(trueData.position_y >= 0));
    trueData.position_x = trueData.position_x(1:row);
    trueData.position_y = trueData.position_y(1:row);
    trueData.velocity_x = trueData.velocity_x(1:row);
    trueData.velocity_y = trueData.velocity_y(1:row);
    n = row;
    t = t(1:row);
end


measData.range = sqrt((trueData.position_x - siteX).^2 + (trueData.position_y - siteY).^2) + ...
    stdRange * randn(n,1);
measData.angle = atan2(trueData.position_y - siteY,trueData.position_x - siteX) + ...
    stdAngle * randn(n,1);
end