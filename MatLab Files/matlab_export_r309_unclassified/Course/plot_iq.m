close all;
clear all;
clc;
addpath('C:\Users\<user>\Documents\MATLAB\matlab2018\Printing');
filename = 'iq_rotator.mp4';
t = 0:0.005:1;
y = exp(1i*2*pi*2.5*t);
k = 5;
pp = PrepForPrint();
figure;
hs1=subplot(1,2,1);
plot(t,[real(y)' imag(y)']);co = get(gca,'ColorOrder');
xlim([-0.1 1.1]);
ylim([-1.1 1.1]);
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
xlabel('Time');
ylabel('Amplitude');

set(gcf,'Color','w');hold on;
hpI1 = plot(t(k),real(y(k)));set(hpI1,'Marker','o','MarkerFaceColor',co(1,:),'MarkerEdgeColor','none');
hpQ1 = plot(t(k),imag(y(k)));set(hpQ1,'Marker','o','MarkerFaceColor',co(2,:),'MarkerEdgeColor','none');
grid on;



hs2=subplot(1,2,2);
a = linspace(-pi,pi,1000);
plot(cos(a),sin(a),'Color',0.7*[1 1 1],'LineWidth',2);hold on;grid on;
set(gca,'DataAspectRatio',[1 1 1],'XLim',[-1.1 1.1],'YLim',[-1.1 1.1]);
hpI2 = plot(real(y(k)),0);set(hpI2,'Marker','o','MarkerFaceColor',co(1,:),'MarkerEdgeColor','none');
hpQ2 = plot(0,imag(y(k)));set(hpQ2,'Marker','o','MarkerFaceColor',co(2,:),'MarkerEdgeColor','none');
hpDot = plot(real(y(k)),imag(y(k)));
set(hpDot,'Marker','o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none');

hLine = line([0 0],[-1.1 1.1],[-1 -1]);
set(hLine,'Parent',hs1,'Color',0.7*[1 1 1],'LineWidth',2);
hLineVert = line([real(y(k)) real(y(k))],[0 imag(y(k))]);
hLineHorz = line([0 real(y(k))],[imag(y(k)) imag(y(k))]);
hLineMag = line([0 real(y(k))],[0 imag(y(k))]);
hLineI = line(0,0);
hLineQ = line(0,0);

set(gcf,'Position',[651.4000  488.2000  876.4000  242.4000]);
%%

set(hLineMag,'Color',[0.5 .5 .5],'LineWidth',2);
set(hLineI,'Color',co(1,:),'LineWidth',2);
set(hLineQ,'Color',co(2,:),'LineWidth',2);
xlabel('Real Axis');
ylabel('Imaginary Axis');
PrepForPrint(gcf,pp);
set([hLineVert hLineHorz],'Color','k','LineStyle',':','LineWidth',1);
 vidobj = VideoWriter(filename,'MPEG-4');
    vidobj.Quality = 15;
    vidobj.FrameRate = 20;
    open(vidobj);
     pause;
for k = 1 : length(t)
    vi = real(y(k));
    vq = imag(y(k));
    
    % td plot
    set(hLine,'XData',[t(k) t(k)]);
    set(hpI1,'XData',t(k),'YData',vi);
    set(hpQ1,'XData',t(k),'YData',vq);
    
    % dots
    set(hpI2,'XData',vi);
    set(hpQ2,'YData',vq);
    set(hpDot,'XData',vi,'YData',vq);
    
    set(hLineVert,'XData',[vi vi],'YData',[0 vq]);
    set(hLineHorz,'XData',[0 vi],'YData',[vq vq]);
    set(hLineMag,'XData',[0 vi],'YData',[0 vq]);
    set(hLineI,'XData',[0 vi],'YData',[0 0]);
    set(hLineQ,'YData',[0 vq],'XData',[0 0]);
    drawnow;
    vidobj.writeVideo(getframe(gcf));

end
close(vidobj);

return;
