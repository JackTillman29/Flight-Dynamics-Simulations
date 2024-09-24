function out =  convertGeodeticWGS84LatLonAltToECEF(lat_deg , lon_deg , alt_m )
WGS84_SEMI_MAJOR = 6378137.0;
WGS84_FIRST_ECCENTRICITY_SQUARED = 0.00669437999014;

cos_lat = cosd(lat_deg);
sin_lat = sind(lat_deg);
cos_lon = cosd(lon_deg);
sin_lon = sind(lon_deg);

%Compute local earth radios based upon WGS-84 datum
Re = WGS84_SEMI_MAJOR ./ sqrt(1.0 - WGS84_FIRST_ECCENTRICITY_SQUARED .* sin_lat .* sin_lat);

%Compute and return ECEF
ecef_x_out_m = (Re + alt_m) .* cos_lat .* cos_lon;
ecef_y_out_m = (Re + alt_m) .* cos_lat .* sin_lon;
ecef_z_out_m = (Re * (1.0 - WGS84_FIRST_ECCENTRICITY_SQUARED) + alt_m) .* sin_lat;

out = [ecef_x_out_m ecef_y_out_m ecef_z_out_m];

end