% Notional aircraft data
x = ReadGPS_Logger('20181022-075351.txt');
ac.LAT84 = x.lat84_deg;
ac.LON84 = x.lon84_deg;
ac.ALT_M = x.alt_m;
ac.HDG   = x.bearing_deg;
ac.ROLL  = 0 * x.bearing_deg; % no data in file
ac.PITCH = 0 * x.bearing_deg; % no data in file

% Notional site data
site.LAT84 = mean(x.lat84_deg);
site.LON84 = mean(x.lon84_deg);
site.ALT_M = 0.0;

% NOTE: NED refers to point on ground beneath a/c in NED
%       SITE refers to point on ground beneath site in NED
%       BODY refers to a/c body frame
% First, get the NED ECEF basis for the ac
ac.B_ECEF2NED = ComputeWGS84_ECEF_Basis(ac.LAT84.',ac.LON84.');

% Next, get the NED ECEF basis for the site
site.B_ECEF2SITE = ComputeWGS84_ECEF_Basis(site.LAT84.',site.LON84.');

% Compute the NED to body matrix for the AC
ac.NED2BODY = RotationMatrixYPR(ac.HDG * pi/180, ac.PITCH * pi/180, ac.ROLL * pi/180);

% This should map the a/c frame to the NED frame of the site
FinalTMX = zeros(3,3,length(ac.LAT84));
for k = 1 : length(ac.LAT84)
    FinalTMX(:,:,k) = site.B_ECEF2SITE * ac.B_ECEF2NED(:,:,k)' * ac.NED2BODY(:,:,k)';
end
ac.yaw_site = zeros(1,length(ac.HDG));
ac.pitch_site = zeros(1,length(ac.PITCH));
ac.roll_site = zeros(1,length(ac.ROLL));
for k = 1 : length(ac.HDG)
    [ac.yaw_site(k),ac.pitch_site(k),ac.roll_site(k)]=Tmx321_to_YPR(FinalTMX(:,:,k)');
end
