close all;
clear;
clc;
figure;
subplot(1,2,1);

 [Xref,Yref,Zref] = sphere(30);
 h_refsphere = surf(Xref,Yref,Zref);
 axis equal;
 set(h_refsphere,'FaceColor','b','EdgeColor','k','EdgeAlpha',0.1,'FaceAlpha',.2);
 
 % array plane
 h_array = patch([-1 -1 1 1],[-1 1 1 -1],[0 0 0 0],'r');
 set(h_array,'FaceAlpha',.5,'FaceColor','k');
 
 h_array_cardinalx = line([-1 1],[0 0]);
 h_array_cardinaly = line([0 0],[-1 1]);
 set([h_array_cardinalx h_array_cardinaly],'Color','k','LineStyle',':');
 
 xlabel('v');ylabel('u');zlabel('w');grid off;set(gca,'XDir','reverse');
 
 az = 30*pi/180;
 el = 0*pi/180;
 [u,v,w] = azel2uvw(az,el);
 h_steerVector = line([0 v],[0 u],[0 w]);
 
 [xcref,ycref]=circle_patch(0,0,0.1,30);
 zc = 0*xcref;
 xc = xcref + v;
 yc = ycref + u;
 
 h_beamw = patch(xc,yc,zc,'r');
 set(h_beamw,'FaceColor','none','EdgeColor','r');
 
%  phi = atan2(xc,yc);
%  thetsq = xc.^2 + yc.^2;
%  thet = -sqrt(thetsq);
 
 zcproj = sqrt(1 - xc.^2 - yc.^2);
 
 
 h_beamw_sphere = patch(xc,yc,zcproj,'r');
 
 % pretty
 set([h_beamw_sphere h_beamw],'FaceColor','r','FaceAlpha',0.5,'EdgeAlpha',0.2);
 
 
 [xaz,yel]=uvw2azel(yc,xc,zcproj);
 subplot(1,2,2);
 h_azel_patch = patch(180/pi*xaz,180/pi*yel,'r');
 axis equal;
 xlim([-90,90]);
 ylim([-90,90]);
 xlabel('Azimuth (deg)');
 ylabel('Elevation (deg)');
 grid on;
 



 for steer_az = -60:3:60
     for steer_el = -30:3:30
         az = steer_az*pi/180;
         el = steer_el*pi/180;
         [us,vs,ws] = azel2uvw(az,el);       
         set(h_steerVector,'XData',[0 vs],'YData',[0 us],'ZData',[0 ws]);
         zc = 0*xcref;
         xc = xcref + vs;
         yc = ycref + us;
         
         set(h_beamw,'XData',xc,'YData',yc);
         zcproj = sqrt(1 - xc.^2 - yc.^2);
         set(h_beamw_sphere,'XData',xc,'YData',yc,'ZData',zcproj);
         
         [xaz,yel]=uvw2azel(yc,xc,zcproj);
         set(h_azel_patch,'XData',180/pi*xaz,'YData',180/pi*yel);
         
         drawnow;
         pause(0.05);
       
     end
 end