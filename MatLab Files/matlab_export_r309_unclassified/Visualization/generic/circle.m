function hLineArray = circle(radius,metersPerPixel,nSegments, center )
angles = 0:(2*pi/nSegments):(2*pi);
angles = [angles 2*pi];
xvals = radius./metersPerPixel * cos(angles) + center(1);
yvals = radius./metersPerPixel * sin(angles) + center(2);
hLineArray = line(xvals,yvals);
end