function out = convertECEF_toTopocentricENU_GeodeticNormal(...
    local_origin_ecef_x, local_origin_ecef_y, local_origin_ecef_z, ...
    geodeticLatDeg,      geodeticLonDeg, ...
    ecef_x_m,  ecef_y_m, ecef_z_m)
    
%Note first set is ECEF origin of local ENU coordinate system
%second is same point in WGS84 lat, lon (no alt)
%third is point of interest in ECEF

sin_phip = sind(geodeticLatDeg);
cos_phip = cosd(geodeticLatDeg);
sin_lambda = sind(geodeticLonDeg);
cos_lambda = cosd(geodeticLonDeg);

ecef_rel_x = (ecef_x_m - local_origin_ecef_x);
ecef_rel_y = (ecef_y_m - local_origin_ecef_y);
ecef_rel_z = (ecef_z_m - local_origin_ecef_z);

%these are from Matlab
east_m  = -sin_lambda .* ecef_rel_x + cos_lambda .* ecef_rel_y;
north_m = -sin_phip .* cos_lambda .* ecef_rel_x + -sin_phip .* sin_lambda .* ecef_rel_y + cos_phip .* ecef_rel_z;
up_m    = cos_phip .* cos_lambda .* ecef_rel_x + cos_phip .* sin_lambda .* ecef_rel_y + sin_phip .* ecef_rel_z;

out = [east_m,north_m,up_m];

end