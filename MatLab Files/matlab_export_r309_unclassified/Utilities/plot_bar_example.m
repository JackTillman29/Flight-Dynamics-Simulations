xpts = [1 2 3 4];
y1   = [1 1.2 1.4 0.8];
y2   = [3 3.5 2.0 6.0];
y3   = [5 7 5.5 8];

figure;
opengl software;
set(gcf,'Renderer','OpenGL');
for k = 1 : length(xpts)
    %PlotBar(xcenter,xwidth,ymin,ymax,varargin)    
    PlotBar(xpts(k),1.0,y1(k),y2(k),'FaceColor','b','FaceAlpha',0.8,'EdgeColor','g');
    PlotBar(xpts(k),1.0,y2(k),y3(k),'FaceColor','r','FaceAlpha',.5);
end

set(gca,'XTick',[1 2 3 4],'XTickLabel',{'a','b','c','d'});