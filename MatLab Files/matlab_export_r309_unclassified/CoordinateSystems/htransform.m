% homogeneous transform
d2r = 180.0 / pi;

% Define Frame A
frA.origin = [0 0 0]';
frA.tmx321 = RotationMatrixYPR(15*d2r,30*d2r,40*d2r);

% Define Frame B
frB.origin = [100 200 300]';
frB.tmx321 = RotationMatrixYPR(-10*d2r,5*d2r,30*d2r);
frB.point = [7 10 20]';

frA.point = (frB.tmx321)' * frB.point + frB.origin;

% now form homogeneous xform
H = zeros(4,4);
H(1:3,1:3) = frB.tmx321';
H(1:3,4)   = frB.origin;
H(4,4)     = 1;

disp('old way');

