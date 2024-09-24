function [ x2,y2 ] = gen_parallelogram_geometry( ...
    angle,length_on_angle,spacing_on_angle, ...
    pattern_offset_along_angle, ...
    length_on_bottom,spacing_on_bottom )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    x = [0:spacing_on_angle:length_on_angle];
    y = 0*x;
    a = angle;
    R = [ cos(a) -sin(a);
          sin(a) cos(a) ];

    rot_el = R*[x;y];
    x = rot_el(1,:);
    y = rot_el(2,:);

    x_bottom = 0:spacing_on_bottom:length_on_bottom;

    nc = length(x_bottom);
    nx = length(x);
    x2 = zeros(length(x)*nc,1);
    y2 = zeros(length(y)*nc,1);
    offset = spacing_on_bottom;
    iPatOffset = 1;
    nPatOffset = length(pattern_offset_along_angle);
    for k = 1 : nc
        %disp(['column ' num2str(k)])
        r_offset = pattern_offset_along_angle(iPatOffset);
        xy_offset = R * [r_offset 0]';
        
        
        x2((1+(k-1)*nx):(k*nx)) = x + (k-1) * offset + xy_offset(1);
        y2((1+(k-1)*nx):(k*nx)) = y + xy_offset(2);
        if(iPatOffset < nPatOffset)
            iPatOffset = iPatOffset + 1;
        else
            iPatOffset = 1;
        end
    end

end

