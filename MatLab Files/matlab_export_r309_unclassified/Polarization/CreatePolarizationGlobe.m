function varargout = CreatePolarizationGlobe(varargin)

%% INPUTS
dataColor = 'b';
if(nargin == 1)
    hax = varargin{1};
    job = [];
    plotall = 1;
elseif(nargin == 2)
    % axes to plot in
    hax = varargin{1};
    % jones vector to plot on globe
    job = varargin{2};
    % select whether or not to plot reference graphics
    plotall = 1;
elseif(nargin == 3)
    % axes to plot in
    hax = varargin{1};
    % jones vector to plot on globe
    job = varargin{2};
    % select whether or not to plot reference graphics
    plotall = varargin{3};
elseif(nargin == 4)
    % axes to plot in
    hax = varargin{1};
    % jones vector to plot on globe
    job = varargin{2};
    % select whether or not to plot reference graphics
    plotall = varargin{3};
    % data color
    dataColor = varargin{4};
end

axes(hax);

set(gca,'XTick',[],'YTick',[],'ZTick',[]);
set(gca,'XColor','w','YColor','w','ZColor','w');

refcolor = 0.7*ones(1,3);
% refcolor = 'b';


%% Plot the Jones Vector on the globe
if(~isempty(job))
r = 1.01;

tmpidx = find(job.delta > pi/2);
job.delta(tmpidx) = job.delta(tmpidx) - pi;
% if(job.delta > pi/2)
%     job.delta = job.delta - pi;
% end

% compv = job.Ex + job.Ey;
% magv = sqrt(real(job.Ex).^2 + real(job.Ey).^2);
% az = atan2(real(job.Ey),real(job.Ex));
% el = atan2(imag(compv),magv)*2;  % mult by 2 to get el to +/- 90deg
% 
% % az = unwrap(az);
% % el = unwrap(el);
% 
% xr = r.*cos(az).*cos(el);
% yr = r.*sin(az).*cos(el);
% zr = r.*sin(el);
% 
% plot3([0 xr],[0 yr],[0 zr],'Color',dataColor,'LineWidth',2)
% text(xr,yr,zr,'Data','VerticalAlignment','top')


% % A = job.A;
% % B = job.B;
% % delta = job.delta;
% % 
% % % compute Stokes Vector
% % S0 = A.^2 + B.^2;
% % S1 = A.^2 - B.^2;
% % S2 = 2.*A.*B.*cos(delta);
% % S3 = 2.*A.*B.*sin(delta);

S0 = 1;
S1 = job.S1;
S2 = job.S2;
S3 = job.S3;

% hs = scatter3(r*S1(:),r*S2(:),r*S3(:),30,'filled','MarkerFaceColor',dataColor);

if(size(dataColor,1) == 1)
    hs = scatter3(r*S1(:),r*S2(:),r*S3(:),30,'filled','MarkerFaceColor',dataColor);
else
    hs = scatter3(r*S1(:),r*S2(:),r*S3(:),30,dataColor,'filled');
end

if(nargout == 1)
    varargout{1} = hs;
end

end

%% PLOT REFERENCE GRAPHICS
if(plotall)
numface = 50;
set(hax,'NextPlot','add');
[xg,yg,zg] = sphere(numface);

faceAlpha = 0.8;
surf(xg,yg,zg,'EdgeColor','none',...
    'FaceColor',0.9*ones(1,3),'FaceAlpha',faceAlpha);

axis equal;
% return

% draw lines of latitude
numpts = 200;
az = linspace(-pi,pi,numpts);
el = linspace(-pi/2,pi/2,10);
r = 1;
for iel = 1:length(el)
xg = r.*cos(az).*cos(el(iel));
yg = r.*sin(az).*cos(el(iel));
zg = r.*sin(el(iel)).*ones(1,numpts);
plot3(xg,yg,zg,'Color',refcolor)
end

% draw lines of longitude
numpts = 200;
el = linspace(-pi/2,pi/2,numpts);
az = linspace(-pi,pi,9)
r = 1;
for iaz = 1:length(az)
xg = r.*cos(az(iaz)).*cos(el);
yg = r.*sin(az(iaz)).*cos(el);
zg = r.*sin(el);
plot3(xg,yg,zg,'Color',refcolor);
end


% draw reference lines
hAlign = 'center'
vAlign = 'middle'
textoffset = 1.1;
refcolor = 0.4*ones(1,3);

% XXXXXX"north" pole (RHC)
% CORRECTION (I think Wikipedia's Stokes Parameters page is wrong)
%  wiki says "RHC = [1 0 0 1]"
%  but math seems to indicate that "RHC = [1 0 0 -1]"
% "south" pole (RHC)
xr = [0 0];
yr = [0 0];
zr = [-1 -1.5];
plot3(xr,yr,zr,'Color',refcolor,'LineWidth',2);
text(xr(2)*textoffset,yr(2)*textoffset,zr(2)*textoffset,'RHC',...
    'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);

% XXXXXX"south" pole (LHC)
% CORRECTION (I think Wikipedia's Stokes Parameters page is wrong)
%  wiki says "RHC = [1 0 0 1]"
%  but math seems to indicate that "RHC = [1 0 0 -1]"
% "north" pole (LHC)
xr = [0 0];
yr = [0 0];
zr = [1 1.5];
plot3(xr,yr,zr,'Color',refcolor,'LineWidth',2);
text(xr(2)*textoffset,yr(2)*textoffset,zr(2)*textoffset,'LHC',...
    'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);

% H-pol
xr = [1 1.5];
yr = [0 0];
zr = [0 0];
plot3(xr,yr,zr,'Color',refcolor,'LineWidth',2);
text(xr(2)*textoffset,yr(2)*textoffset,zr(2)*textoffset,'H',...
    'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);

% V-pol
xr = -1*[1 1.5];
yr = [0 0];
zr = [0 0];
plot3(xr,yr,zr,'Color',refcolor,'LineWidth',2);
text(xr(2)*textoffset,yr(2)*textoffset,zr(2)*textoffset,'V',...
    'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);


% slant 45deg
xr = [0 0];
yr = [1 1.5];
zr = [0 0];
plot3(xr,yr,zr,'Color',refcolor,'LineWidth',2);
text(xr(2)*textoffset,yr(2)*textoffset,zr(2)*textoffset,'S.45',...
    'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);
% slant 135deg
xr = [0 0];
yr = -1*[1 1.5];
zr = [0 0];
plot3(xr,yr,zr,'Color',refcolor,'LineWidth',2);
text(xr(2)*textoffset,yr(2)*textoffset,zr(2)*textoffset,'S.135',...
    'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);

view([135 30])


end

return


%% H
scale = 0.1;
xc = scale*[-1 1];
yc = [0 0];
xl = 0;
yl = 0;
plot(xc+xl,yc+yl,'Color',refcolor)

%% V
xc = scale*[0 0];
yc = scale*[-1 1];
xl = -1;
yl = 0;
plot(xc+xl,yc+yl,'Color',refcolor)

% if mercadian proj
xc = scale*[0 0];
yc = scale*[-1 1];
xl = 1;
yl = 0;
plot(xc+xl,yc+yl,'Color',refcolor)

%% Slant +45deg
xc = scale*[-1 1];
yc = scale*[-1 1];
xl = 0.5;
yl = 0;
plot(xc+xl,yc+yl,'Color',refcolor)

%% Slant -45deg
xc = scale*[-1 1];
yc = scale*[1 -1];
xl = -0.5;
yl = 0;
plot(xc+xl,yc+yl,'Color',refcolor)


%% CIRCLE TEMPLATE

numpts = 100;
ileft = 1;
itop = 25;
iright = 50;
ibot = 75;

th = linspace(-pi,pi,numpts);
scale = 0.1;
xc = scale*cos(th);
yc = scale*sin(th);

%% LHC ref
for xl = [-1:0.5:1]
yl = 1
plot(xc+xl,yc+yl,'Color',refcolor)
scatter(xc(1)+xl,yc(1)+yl,25,'filled',...
    'Marker','^','MarkerFaceColor',refcolor)
end

%% RHC ref
for xl = [-1:0.5:1]
yl = -1
plot(xc+xl,yc+yl,'Color',refcolor)
scatter(xc(1)+xl,yc(1)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)
end

%% LH elliptical
xl = 0
yl = 0.5
xe = 1.0
ye = 0.5
plot(xe*xc+xl,ye*yc+yl,'Color',refcolor)
scatter(xe*xc(1)+xl,ye*yc(1)+yl,25,'filled',...
    'Marker','^','MarkerFaceColor',refcolor)

%% RH elliptical
xl = 0
yl = -0.5
xe = 1.0
ye = 0.5
plot(xe*xc+xl,ye*yc+yl,'Color',refcolor)
scatter(xe*xc(1)+xl,ye*yc(1)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)

%% slant +45 LH elliptical
xl = 0.5
yl = 0.5
xe = 1.0
ye = 0.5
rotang = 45;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
scatter(xce(ileft)+xl,yce(ileft)+yl,25,'filled',...
    'Marker','^','MarkerFaceColor',refcolor)

%% slant +90 LH elliptical
xl = 1.0
yl = 0.5
xe = 1.0
ye = 0.5
rotang =90;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
hc = scatter(xce(ibot)+xl,yce(ibot)+yl,25,'filled',...
    'Marker','^','MarkerFaceColor',refcolor)

%% slant -45 LH elliptical
xl = -0.5
yl = 0.5
xe = 1.0
ye = 0.5
rotang = -45;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
scatter(xce(ileft)+xl,yce(ileft)+yl,25,'filled',...
    'Marker','^','MarkerFaceColor',refcolor)

%% slant -90 LH elliptical
xl = -1.0
yl = 0.5
xe = 1.0
ye = 0.5
rotang = -90;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
scatter(xce(itop)+xl,yce(itop)+yl,25,'filled',...
    'Marker','^','MarkerFaceColor',refcolor)











%% slant +45 RH elliptical
xl = 0.5
yl = -0.5
xe = 1.0
ye = 0.5
rotang = 45;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
scatter(xce(ileft)+xl,yce(ileft)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)

%% slant +90 RH elliptical
xl = 1.0
yl = -0.5
xe = 1.0
ye = 0.5
rotang =90;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
hc = scatter(xce(ibot)+xl,yce(ibot)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)

%% slant -45 RH elliptical
xl = -0.5
yl = -0.5
xe = 1.0
ye = 0.5
rotang = -45;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
scatter(xce(ileft)+xl,yce(ileft)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)

%% slant -90 RH elliptical
xl = -1.0
yl = -0.5
xe = 1.0
ye = 0.5
rotang = -90;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
hc = scatter(xce(itop)+xl,yce(itop)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)


%% Assign tick labels

% ellipticity angles in deg (normalized plotting above)
%  vary between +/- 45deg
set(gca,'XTick',[-1 -0.5 0 0.5 1])
set(gca,'XTickLabel',{'45^o','22.5^o','0^o','-22.5^o','-45^o'})

% rotation angles in deg (normalized plotting above)
%   vary between +/- 90deg
set(gca,'YTick',[-1 -0.5 0 0.5 1])
set(gca,'YTickLabel',{'90^o','45^o','0^o','-45^o','-90^o'})


end