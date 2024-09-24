function CreatePolarizationMap(varargin)
% function CreatePolarizationMap(hax,jones_obj,plotall,dataColor)

% figure
% hax = axes('NextPlot','add');

%% INPUTS
dataColor = 'b';
if(nargin == 1)
    hax = varargin{1};
    jones_obj = [];
    plotall = 1;
elseif(nargin == 2)
    % axes to plot in
    hax = varargin{1};
    % jones vector to plot on globe
    jones_obj = varargin{2};
    % select whether or not to plot reference graphics
    plotall = 1;
elseif(nargin == 3)
    % axes to plot in
    hax = varargin{1};
    % jones vector to plot on globe
    jones_obj = varargin{2};
    % select whether or not to plot reference graphics
    plotall = varargin{3};
elseif(nargin == 4)
    % axes to plot in
    hax = varargin{1};
    % jones vector to plot on globe
    jones_obj = varargin{2};
    % select whether or not to plot reference graphics
    plotall = varargin{3};
    % data color
    dataColor = varargin{4};
end


axes(hax)

refcolor = 0.7*[1 1 1];
% refcolor = 'b';

if(plotall)
%% H
scale = 0.1;
xc = scale*[-1 1];
yc = [0 0];
xl = 0;
yl = 0;
set(hax,'NextPlot','add','Box','On')
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
rotang = -45;
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
rotang = -90;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
hc = scatter(xce(itop)+xl,yce(itop)+yl,25,'filled',...
    'Marker','v','MarkerFaceColor',refcolor)

%% slant -45 RH elliptical
xl = -0.5
yl = -0.5
xe = 1.0
ye = 0.5
rotang = 45;
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
rotang = 90;
    xce = xe.*xc.*cosd(rotang) - ye*yc.*sind(rotang);
    yce = xe.*xc.*sind(rotang) + ye*yc.*cosd(rotang);
plot(xce+xl,yce+yl,'Color',refcolor)
hc = scatter(xce(ibot)+xl,yce(ibot)+yl,25,'filled',...
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

% hp_empty = plot([0 0],[0 0],'w','linewidth',0.0001);
% legend(hp_empty,'Propagation Direction Out of the Page')
title('Propagation Direction Out of the Page')
end

%% Plot data
% % compv = jones_obj.Ex + jones_obj.Ey;
% % magv = sqrt(real(jones_obj.Ex).^2 + real(jones_obj.Ey).^2);
% % az = atan2(real(jones_obj.Ey),real(jones_obj.Ex));
% % el = atan2(imag(compv),magv)*2;  % mult by 2 to get el to +/- 90deg
% % 
% % az = 2*az/pi;
% % el = 2*el/pi;
% % 
% % if(az > 1 & az < 2)
% %     az = mod(az,1)-1;
% % elseif(az > -1 & az < 0)
% %     az = mod(az,1)-1;
% % else
% %     az = mod(az,1);
% % end
% % 
% % % az = mod(az,1);
% % el = mod(el,1);

if(~isempty(jones_obj))
% A = jones_obj.A;
% B = jones_obj.B;
% delta = jones_obj.delta;
% 
% % compute Stokes Vector
% S0 = A.^2 + B.^2;
% S1 = A.^2 - B.^2;
% S2 = 2.*A.*B.*cos(delta);
% S3 = 2.*A.*B.*sin(delta);
% 
% if(abs(S1) < 1e-6)
%     S1 = 0;
% end
% 
% if(abs(S2) < 1e-6)
%     S2 = 0;
% end
% 
% if(abs(S3) < 1e-6)
%     S3 = 0;
% end

[az,el,r] = cart2sph(jones_obj.S1(:),jones_obj.S2(:),jones_obj.S3(:));
az = az/pi * 1;
el = el/pi * 2;

% % x = jones_obj.alpha/(pi/4)/2;
% % 
% % if(jones_obj.delta > pi/2)
% %     y = 0
% % else
% %     y = jones_obj.delta/pi;
% % end
% % scatter(x,y,30,'filled')

if(size(dataColor,1) == 1)
    scatter(az,el,30,'filled','MarkerFaceColor',dataColor)
else
    scatter(az,el,30,dataColor,'filled')
end

xlabel('Ellipticity Angle [deg]')
ylabel('Rotation Angle [deg]')

xlim(1.2*[-1 1])
ylim(1.2*[-1 1])
end

end