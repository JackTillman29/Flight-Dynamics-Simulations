function [B,varargout] = ComputeWGS84_ECEF_Basis(lat84,lon84)
% function [B,varargout] = ComputeWGS84_ECEF_Basis(lat84,lon84)
% Inputs: Latitude (deg) in WGS84 Datum (Nx1)
%         Longitude (deg) in WGS84 Datum (Nx1)
% Outputs: 
%      [ unit vector of increasing North in ECEF ] 
%  B = [ unit vector of increasing East  in ECEF ] 
%      [ unit vector of increasing Down  in ECEF ] 
%  varargout{1} = uN (3x1)
%  varargout{2} = uE (3x1)
%  varargout{3} = uD (3x1)
%  
%  Note: B forms the Local ECEF->NED transformation matrix, that is:
%        x_ned = B * x_ecef
% K. Sawmiller 2018
uN = [-sind(lat84).*cosd(lon84) -sind(lat84).*sind(lon84) cosd(lat84)];
uD = -[cosd(lat84).*cosd(lon84) cosd(lat84).*sind(lon84) sind(lat84)];
uE = [-sind(lon84) cosd(lon84) 0*lon84];
B = zeros(3,3,length(lat84));
%B = [uN;uE;uD];
for k = 1 : length(lat84)
    B(:,:,k) = [uN(k,:);uE(k,:);uD(k,:)];
end
if(nargout > 1)
    varargout{1}=uN';
    varargout{2}=uE';
    varargout{3}=uD';
end

end