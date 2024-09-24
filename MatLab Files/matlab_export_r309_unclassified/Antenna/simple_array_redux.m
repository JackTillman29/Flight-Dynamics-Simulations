close all;
clear;
clc;



% Create planar array element positions
NeY = 20; % # of elements
NeZ = 20;
Ne = NeY * NeZ;


[yG,zG] = meshgrid( ...
    linspace(1,NeY,NeY), ...
    linspace(1,NeZ,NeZ) );

y = yG(:);
z = zG(:);
x = 0 * y;

y = y - mean(y);
z = z - mean(z);

% 150MHz, wavelength = 2m, lambda/2 = 1m
lambda = 3e8 / 150e6;
azSteer = 0*pi/180;
elSteer = 0;
[uxs,uys,uzs]=sph2cart(azSteer,elSteer,1);

p = [x y z];

% compute far-field evaluation locations (# of locations = M)
nAz = 200;
nEl = 200;
azPts = linspace(-pi/2,pi/2,nAz);
elPts = linspace(-pi/2,pi/2,nEl);
[azG,elG]=meshgrid( ...
    azPts, ...
    elPts ...
    );

% convert to unit vector
[ux,uy,uz]=sph2cart(azG,elG,1);

R = [ ux(:)'
      uy(:)'
      uz(:)' ];
  
rs = [uxs uys uzs]';
  
phs_steer = -2*pi / lambda * p * rs;

% build up random phase 
phsnz_mag = 5;
phsnz = phsnz_mag*pi/180*randn(Ne,1);

% Compute phase map [N x 3] [3 x M]
P = 2*pi / lambda * p * R ;

% Apply weighing to complex exponential and sum
% Scale based on # of elements
wa = ones(Ne,1);
ws = exp(1i*(phs_steer + phsnz));

w = wa .* ws;
PW = w.' * exp(1i*P) ./ sqrt(Ne);
Res = reshape(PW,nEl,nAz);

%% Plotting 
hi=imagesc(elPts * 180/pi, azPts * 180/pi, 20*log10(abs(Res)));
set(gca,'YDir','normal');
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
colorbar;
caxis([-40 40]);

set(hi,'ButtonDownFcn','patternBtnDown');


[px,py,pz] = sphere(50);

figure;
sEarth = surface(py, px ,(pz));  
sEarth.FaceColor = 'texturemap';        % set color to texture mapping
sEarth.EdgeColor = 'none';              % remove surface edge color
sEarth.CData = [get(hi,'CData') 0*get(hi,'CData')-99];

set(gca,'DataAspectRatio',[1 1 1],'Color','k','XColor','w','YColor','w','ZColor','w');
ylabel('U Space');
zlabel('V Space');
hc=colorbar;
set(hc,'Color','w')
caxis([-40 40]);
set(gcf,'Color','k');

