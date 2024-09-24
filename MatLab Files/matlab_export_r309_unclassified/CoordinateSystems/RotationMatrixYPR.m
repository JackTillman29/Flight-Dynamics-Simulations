function R = RotationMatrixYPR(y,p,r,angmode)
%        Syntax: R = RotationMatrixYPR(y,p,r,angmode)
% Input:
%        Yaw   (rad/deg) from A (inertial) to B (body)
%        Pitch (rad/deg) from A (inertial) to B (body)
%        Roll  (rad/deg) from A (inertial) to B (body)
%        **NOTE**: Code assumes a Euler 3,2,1 (Yaw,Pitch,Roll) sequence/
%
%        angmode: 'rad' = radians, 'deg' = degrees
%
% Output:
%        3x3 Rotation matrix             Projection of basis A onto B
%        [--- B basis #1 in A ----]      [A1.B1 A2.B1 A3.B1]
%        [--- B basis #2 in A ----]  or  [A1.B2 A2.B2 A3.B2]
%        [--- B basis #3 in A ----]      [A1.B3 A2.B3 A3.B3]
%
%                                 equivalently.....
%
%        [         |                     |                  |
%        [  A basis #1 in B      A basis #2 in B    A basis #3 in B ]
%        [         |                     |                  |
%
% Linear operator projects vectors defined in basis "A" onto basis "B"
% Transposed operator projects vectors defined in basis "B" onto basis "A"
% K. Sawmiller, 11/2018

if(~exist('angmode','var'))
    angmode = 'rad';
end
switch(lower(angmode))
    case 'rad'
        sf = 1.0;
    case 'deg'
        sf = pi/180;
    otherwise
        error(['Unknown option: ' angmode]);
end

cy = cos(sf*y);
cp = cos(sf*p);
cr = cos(sf*r);
sy = sin(sf*y);
sp = sin(sf*p);
sr = sin(sf*r);
R = zeros(3,3,length(y));

R(1,1,:) = cp.*cy;
R(1,2,:) = cp.*sy;
R(1,3,:) = -sp;

R(2,1,:) = sr.*sp.*cy - cr.*sy;
R(2,2,:) = cr.*cy + sr.*sp.*sy;
R(2,3,:) = sr.*cp;

R(3,1,:) = sr.*sy + cr.*sp.*cy;
R(3,2,:) = cr.*sp.*sy - sr.*cy;
R(3,3,:) = cr.*cp;


end