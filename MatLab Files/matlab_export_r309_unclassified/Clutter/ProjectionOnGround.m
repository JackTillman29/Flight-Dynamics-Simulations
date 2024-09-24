function [vbwGndRng,hbwGndRng,varargout]= ProjectionOnGround(slantRngToGnd,angOffHorz_rad,beamwidth_rad)

d2r = pi/180;
r2d = 1/d2r;

alt = slantRngToGnd .* cos(angOffHorz_rad);
gndRng = slantRngToGnd .* sin(angOffHorz_rad);

R = slantRngToGnd;
th = angOffHorz_rad;
beta = beamwidth_rad/2;

%==========================================================================
% GET PROJECTED VERTICAL BEAMWIDTH RANGE EXTENT ON GROUND
h = R * tan(beta);

% inner leg of beamwidth projection onto ground
a  = 90*d2r - th; % angle OPPOSITE h in non-right-angle triangle
b1 = 90*d2r - beta;
H1 = 180*d2r - a - b1;
vbwGndRng1 = h .* sin(b1) ./ sin(H1);

% outer leg
b2 = 180*d2r - b1;
H2 = 180*d2r - a - b2;
vbwGndRng2 = h .* sin(b2) ./ sin(H2);

% phi   = b1*r2d
% gamma = H1*r2d
% xi    = b2*r2d
% kappa = H2*r2d

vbwGndRng = vbwGndRng1 + vbwGndRng2;

%==========================================================================
% GET PROJECTED HORIZONTAL BEAMWIDTH RANGE EXTENT ON GROUND

hbwGndRng1 = h;
hbwGndRng2 = h;
hbwGndRng  = hbwGndRng1 + hbwGndRng2;


if(nargout==3)
    data.vbwRng1 = vbwGndRng1;
    data.vbwRng2 = vbwGndRng2;
    data.vbwRng  = vbwGndRng;
    
    data.hbwRng1 = hbwGndRng1;
    data.hbwRng2 = hbwGndRng2;
    data.hbwRng  = hbwGndRng;
    
    varargout{1} = data;
end

end