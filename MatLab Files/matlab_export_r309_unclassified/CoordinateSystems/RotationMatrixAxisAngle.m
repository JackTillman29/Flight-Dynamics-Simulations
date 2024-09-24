function R = RotationMatrixAxisAngle(a,angRad)
    % Define vector a as a 3x1 column vector, unit magnitude
    % R = cos(angRad) * I + (1 - cos(angRad)) * [a * a'] + sin(angRad) * [0 a(3) -a(2);-a(3) 0 a(1);a(2) -a(1) 0]
    R = cos(angRad)*eye(3) + (1 - cos(angRad)) * (a * a') + sin(angRad) * [0 a(3) -a(2);-a(3) 0 a(1);a(2) -a(1) 0];
end