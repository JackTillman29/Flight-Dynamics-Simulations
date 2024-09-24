function draw_polarization_in_antenna_frame(txpos,Ex,Ey,txframe_x,txframe_y,txframe_z,varargin)

if(length(varargin) > 0)
    hax = varargin{1};
end

R = [reshape(txframe_x,[3 1]), ...
     reshape(txframe_y,[3 1]), ...
     reshape(txframe_z,[3 1])];

nPerWave = 100;
rPerN    = 2*pi./nPerWave;
nWaves = 5;
kvec = 0:(nPerWave*nWaves);

% COMPUTE JONES VECTOR
Eeff = sqrt(abs(Ex).^2 + abs(Ey).^2);
   Ex = Ex./Eeff;
   Ey = Ey./Eeff;
   A = abs(Ex)./Eeff;
   B = abs(Ey)./Eeff;
   deltay = angle(Ey);
   deltax = angle(Ex);
   delta  = deltay - deltax;

nPerWave = 100;
rPerN    = 2*pi./nPerWave;
nWaves = 5;
kvec = 0:(nPerWave*nWaves);

%A represents x,z(k) plane
x = A * cos(rPerN * kvec);

%B represents y,z(k) plane rel x
y = B * cos(rPerN * kvec + delta);

% quiver x to show DIRECTION of rotation
    karrow = nPerWave/nWaves
    karrowp1 = round(0.1*nPerWave);
    karrowp1 = 1;
    qt = kvec(karrow);
    qx = x(karrow);
    qy = y(karrow);

    dt = kvec(karrow+karrowp1) - kvec(karrow);
    dx = x(karrow+karrowp1) - x(karrow);
    dy = y(karrow+karrowp1) - y(karrow);

    ang = 30;
    RT1 = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
    RT2 = [cosd(-ang) -sind(-ang); sind(-ang) cosd(-ang)];
    rotpt1 = RT1*[-dx -dy].';
    rotpt2 = RT2*[-dx -dy].';
    p1 = qx + 2*[0 rotpt1(1)];
    p2 = qy + 2*[0 rotpt1(2)];
    p3 = qx + 2*[0 rotpt2(1)];
    p4 = qy + 2*[0 rotpt2(2)];
    
    arpt1 = [[0 0]; p1; p2];
    arpt1 = R * arpt1;
    arpt2 = [[0 0]; p3; p4];
    arpt2 = R * arpt2;
    
    p1 = R*[0; p1']; % convert to 3D vectors in the YZ plane and rotate
    p2 = R*[0; p2'];
    p3 = R*[0; p3'];
    p4 = R*[0; p4'];
    
polcurve = [x.*0; x; y];
% figure; plot3(x.*0,x,y); hold on;
% hp1 = plot3([0 0],p1,p2,'k','linewidth',2);
% hp2 = plot3([0 0],p3,p4,'k','linewidth',2);
% return

polcurve_rot = R * polcurve;
x = polcurve_rot(1,:);
y = polcurve_rot(2,:);
z = polcurve_rot(3,:);


if(exist('hax'))
%     plot3(hax, ...
%           txpos(1) + [0 txframe_x(1)], ...
%           txpos(2) + [0 txframe_x(2)], ...
%           txpos(3) + [0 txframe_x(3)], ...
%           'r')
%     hold on; axis equal
%     plot3(hax, ...
%           txpos(1) + [0 txframe_y(1)], ...
%           txpos(2) + [0 txframe_y(2)], ...
%           txpos(3) + [0 txframe_y(3)], ...
%           'g')
%     plot3(hax, ...
%           txpos(1) + [0 txframe_z(1)], ...
%           txpos(2) + [0 txframe_z(2)], ...
%           txpos(3) + [0 txframe_z(3)], ...
%           'b')
    plot3(hax,txpos(1)+x,         txpos(2)+y,         txpos(3)+z,'k')
    plot3(hax,txpos(1)+arpt1(1,:),txpos(2)+arpt1(2,:),txpos(3)+arpt1(3,:),'k','linewidth',2);
    plot3(hax,txpos(1)+arpt2(1,:),txpos(2)+arpt2(2,:),txpos(3)+arpt2(3,:),'k','linewidth',2);
else
%     plot3(txpos(1) + [0 txframe_x(1)], ...
%           txpos(2) + [0 txframe_x(2)], ...
%           txpos(3) + [0 txframe_x(3)], ...
%           'r')
%     hold on; axis equal
%     plot3(txpos(1) + [0 txframe_y(1)], ...
%           txpos(2) + [0 txframe_y(2)], ...
%           txpos(3) + [0 txframe_y(3)], ...
%           'g')
%     plot3(txpos(1) + [0 txframe_z(1)], ...
%           txpos(2) + [0 txframe_z(2)], ...
%           txpos(3) + [0 txframe_z(3)], ...
%           'b')
%     plot3(x,y,z,'k')
    plot3(txpos(1)+arpt1(1,:),txpos(2)+arpt1(2,:),txpos(3)+arpt1(3,:),'k','linewidth',2);
    plot3(txpos(1)+arpt2(1,:),txpos(2)+arpt2(2,:),txpos(3)+arpt2(3,:),'k','linewidth',2);
end

end
