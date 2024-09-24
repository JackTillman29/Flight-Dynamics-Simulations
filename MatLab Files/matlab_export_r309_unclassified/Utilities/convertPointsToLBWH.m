function LBWH = convertPointsToLBWH(pts)
% converts a set of points to the largest surrounding box
% pts: Nx2
% LBWH: 1x4 (left, bottom, width, height)

LBWH(1) = min(pts(:,1));
LBWH(2) = min(pts(:,2));

dist_x = pts(:,1) - pts(:,1).';
dist_y = pts(:,2) - pts(:,2).';

LBWH(3) = max(dist_x(:));
LBWH(4) = max(dist_y(:));

end