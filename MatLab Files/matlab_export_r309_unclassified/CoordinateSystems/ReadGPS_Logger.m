function y = ReadGPS_Logger(infile)
% Author: Keith Sawmiller, 2018
% Purpose: Reads in Android "GPS Logger" data file format. Returns a
% structure with the following fields:
% Time (s)
% Latitude WGS84 (deg)
% Longitude WGS84 (deg)
% Altitude above WGS84 ellipsoid (m)
% Accuracy (m)
% Speed (m/s)
% Bearing (deg)
% # Satellites in view
% # Satellites used for solution
    x=importdata(infile);
    N = length(x.textdata);
    y.t = zeros(1,N-1);
    y.lat84_deg = y.t;
    y.lon84_deg = y.t;
    y.accur_m = y.t;
    y.alt_m = y.t;
    
    y.speed_mps = y.t;
    y.bearing_deg = y.t;
    y.sat_used = y.t;
    y.sat_inview = y.t;
    
    for k = 2 : N
        hms=sscanf(strrep(sscanf(x.textdata{k,2},'%*s %s'),':',' '),'%d')';
        y.t(k-1) = hms(1)*3600+hms(2)*60+hms(3);
    end
    y.lat84_deg(1) = str2num(x.textdata{2,3});
    y.lon84_deg(1) = str2num(x.textdata{2,4});
    y.accur_m(1) = str2num(x.textdata{2,5});
    y.alt_m(1) = str2num(x.textdata{2,6});
    y.speed_mps(1) = str2num(x.textdata{2,8});
    y.bearing_deg(1) = str2num(x.textdata{2,9});
    y.sat_used(1) = str2num(x.textdata{2,10});
    y.sat_inview(1) = str2num(x.textdata{2,11});
    
    y.lat84_deg(2:end) = x.data(:,1);
    y.lon84_deg(2:end) = x.data(:,2);
    y.accur_m(2:end) = x.data(:,3);
    y.alt_m(2:end) = x.data(:,4);
    y.speed_mps(2:end) = x.data(:,6);
    y.bearing_deg(2:end) = x.data(:,7);
    y.sat_used(2:end) = x.data(:,8);
    y.sat_inview(2:end) = x.data(:,9);
end