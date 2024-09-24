function [x,y] = ellipse_patch(centerx,centery,a,b,segments)
%function [x,y] = circle_patch(centerx,centery,a,b,segments)

    % 
    xTop = -a:(2*a/segments):a;
    
    y2 = b^2 * ( 1 - (xTop./a).^2 );
    
    yTop = sqrt(y2);
    
    xBottom = xTop((end-1):-1:2);
    yBottom = -yTop((end-1):-1:2);
    
    x = [xTop xBottom] + centerx;
    y = [yTop yBottom] + centery;

    
end