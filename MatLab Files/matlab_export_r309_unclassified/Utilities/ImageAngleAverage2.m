function [rng_avg_accum,rng_std_accum] = ImageAngleAverage2(xgrid,ygrid,cdata,val,origin,vertang,angwidth,plotflag)
% Purpose of this function is to compute an average aspect value based upon
% the inputs provided
%
% Things we need: 
% 1 - Image CData,xgrid,ygird
% 2 - Value of interest
% 3 - 0,0 point
% 4 - Angle from vertical
% 5 - Aspect width from (4) in degrees
%
% Returns average value

if(~exist('plotflag'))
    plotflag = 0;
end

[XGRID,YGRID] = meshgrid(xgrid,ygrid);
RGRID = sqrt(XGRID.^2 + YGRID.^2);
vert_angle_grid = atan2(XGRID,YGRID);

ang_test_min = vertang-angwidth/2.0;
ang_test_max = vertang+angwidth/2.0;

% argument check. if both angles are > 180, adjust to -angles
if(ang_test_min > pi)
    ang_test_min = ang_test_min - 2*pi;
    ang_test_max = ang_test_max - 2*pi;
end

cdata_up   = cdata([2:end 1],:);
cdata_down = cdata([end 1:(end-1)],:);
cdata_left = cdata(:,[2:end 1]);
cdata_right = cdata(:,[end 1:(end-1)]);

overlay_image = cdata_right + cdata_left + cdata_up + cdata_down;
overlay_image(cdata == 1) = 0;
overlay_image(overlay_image ~= 0) = 1;


% New Stuff
xp = 1000*1852*[0 sin(ang_test_min) sin(ang_test_max)];
yp = 1000*1852*[0 cos(ang_test_min) cos(ang_test_max)];
INPOLY = (overlay_image & inpolygon(XGRID,YGRID,xp,yp));
rTest = RGRID(INPOLY == 1);

rng_avg_accum = mean(rTest);
rng_std_accum = std(rTest);

% debug region
if(nargout == 0 || plotflag == 1)
    figure;
    imagesc(xgrid,ygrid,cdata);
    set(gca,'YDir','normal');
    hold on;
    plot(origin(1),origin(2),'wo');
    axis equal;
    
    n  = 100;
    xp = zeros(1,100);
    yp = zeros(1,100);
    xpplus  = xp;
    xpminus = xp;
    ypplus  = yp;
    ypminus = yp;
    
    astep = angwidth / (n-1);
    for k = 1 : n
        xp(k) = rng_avg_accum * sin(ang_test_min + (k-1) * astep);
        xpplus(k)  = (rng_avg_accum + rng_std_accum)* sin(ang_test_min + (k-1) * astep);
        xpminus(k) = (rng_avg_accum - rng_std_accum)* sin(ang_test_min + (k-1) * astep);
        yp(k) = rng_avg_accum * cos(ang_test_min + (k-1) * astep);
        ypplus(k)  = (rng_avg_accum + rng_std_accum)* cos(ang_test_min + (k-1) * astep);
        ypminus(k) = (rng_avg_accum - rng_std_accum)* cos(ang_test_min + (k-1) * astep);
    end
    
    h1 = line(xp,yp);
    set(h1,'Color','w');
    h2 = line(xpplus,ypplus);
    h3 = line(xpminus,ypminus);
    set([h2 h3],'Color','w','LineStyle',':');

end

end