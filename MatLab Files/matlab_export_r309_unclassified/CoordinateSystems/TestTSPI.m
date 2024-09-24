% Notional aircraft data
ac.LAT84 = 1;
ac.LON84 = 0;
ac.ALT_M = 9144;
ac.HDG   = 90.0;
ac.ROLL  = 0.0;
ac.PITCH = 0.0;

% Notional site data
site.LAT84 = 0;
site.LON84 = 0.0;
site.ALT_M = 0.0;

% NOTE: NED refers to point on ground beneath a/c in NED
%       SITE refers to point on ground beneath site in NED
%       BODY refers to a/c body frame
% First, get the NED ECEF basis for the ac
ac.B_ECEF2NED = ComputeWGS84_ECEF_Basis(ac.LAT84,ac.LON84);

% Next, get the NED ECEF basis for the site
site.B_ECEF2SITE = ComputeWGS84_ECEF_Basis(site.LAT84,site.LON84);

% Compute the NED to body matrix for the AC
ac.NED2BODY = RotationMatrixYPR(ac.HDG * pi/180, ac.PITCH * pi/180, ac.ROLL * pi/180);

% This should map the a/c frame to the NED frame of the site
FinalTMX = site.B_ECEF2SITE * ac.B_ECEF2NED' * ac.NED2BODY';

[y,p,r]=Tmx321_to_YPR(FinalTMX');
disp(['Yaw  : ' num2str(180/pi*y)]);
disp(['Pitch: ' num2str(180/pi*p)]);
disp(['Roll : ' num2str(180/pi*r)]);