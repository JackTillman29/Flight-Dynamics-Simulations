function [x,y] = circle_patch(centerx,centery,radius,segments)
%function [x,y] = circle_patch(centerx,centery,radius,segments)

    % build circle
    seg_ang = 2*pi / segments;
    xyang = 0:seg_ang:(2*pi-seg_ang);
    
    xp = radius * cos(xyang);
    yp = radius * sin(xyang);
    
    x = xp + centerx;
    y = yp + centery;

    
end