function R = RotationMatrixAxisReflection(a)
    % Assumes that the rotation about an axis is 180 degrees, typically
    % used for ray trace applications
    %
    % cos(pi) = -1
    % sin(pi) =  0
    % R = I + (2) * [a' * a]
    R = -eye(3) + 2 * (a * a');
end