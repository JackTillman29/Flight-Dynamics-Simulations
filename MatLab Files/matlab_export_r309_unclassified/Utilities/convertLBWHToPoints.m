function pts = convertLBWHToPoints(LBWH)
% converts a box defined by LBWH (Left,Bottom,Width,Height) into 4 points
% LBWH:  length-4
% pts:   size(4,2)

pts(1,1:2) = [LBWH(1) LBWH(2)];
pts(2,1:2) = [(LBWH(1)+LBWH(3)) LBWH(2)];
pts(3,1:2) = [(LBWH(1)+LBWH(3)) (LBWH(2)+LBWH(4))];
pts(4,1:2) = [LBWH(1) (LBWH(2)+LBWH(4))];


end