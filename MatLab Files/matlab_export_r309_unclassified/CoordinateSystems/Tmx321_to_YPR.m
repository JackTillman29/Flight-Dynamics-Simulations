function [yaw,pitch,roll] = Tmx321_to_YPR(R,angmode)
%          syntax: [yaw,pitch,roll] = Tmx321_to_YPR(R,angmode)
% Input:
%        3x3 Rotation matrix              projection of basis A onto B
%        [--- A basis #1 in B ----]      [A1.B1 A1.B2 A1.B3]
%        [--- A basis #2 in B ----]  or  [A2.B1 A2.B2 A2.B3]
%        [--- A basis #3 in B ----]      [A3.B1 A3.B2 A3.B3]
% Output: 
%        Yaw   (rad/deg) from A (inertial) to B (body)
%        Pitch (rad/deg) from A (inertial) to B (body)
%        Roll  (rad/deg) from A (inertial) to B (body)
%
%        angmode: 'rad' = radians, 'deg' = degrees
%
% Input linear operator projects vectors defined in basis "A" onto basis "B"
% Or simply, it transforms inertial (A) vectors into the body (B) frame
% The output angles are the equivalent Euler angles to achieve the same
% effect
% K. Sawmiller, 11/2018

if(~exist('angmode','var'))
    angmode = 'rad';
end
switch(lower(angmode))
    case 'rad'
        sf = 1.0;
    case 'deg'
        sf = 180/pi;
    otherwise
        error(['Unknown option: ' angmode]);
end

yaw   = atan2(R(1,2),R(1,1));
pitch = asin(-R(1,3));
roll = atan2(R(2,3),R(3,3));

yaw = sf*yaw;
pitch = sf*pitch;
roll = sf*roll;

end

