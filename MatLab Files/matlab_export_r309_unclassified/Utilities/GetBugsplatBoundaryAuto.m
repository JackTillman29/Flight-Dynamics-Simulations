function [dr,cr] = GetBugsplatBoundaryAuto(hax,units)

tmpfig = figure;
hax = copyobj(hax,tmpfig);

% small
% set(gcf,'position',[561  279  584  520]);
% smaller
% set(gcf,'position',[948  382  403  369]);
set(gcf,'position',[948  382  315  292]);

% set colormap to be JET
set(gcf,'Colormap',jet);
% set background color to white
set(gcf,'Color','w')

% undo all axes formatting:
grid off; box off
colorbar off
set(gca,'XTick',[])
set(gca,'YTick',[])

xlabel(hax,'')
ylabel(hax,'')
title(hax,'')
% get the axes position in pixels
set(gca,'units','pixels')
axes_pos = round(get(gca,'position'));

colormap([0 0 0])
% get bitmap
framedat = getframe(gcf);
plotdat  = framedat.cdata(:,:,1);

xlims = get(gca,'XLim')/units;
ylims = get(gca,'YLim')/units;

minx = xlims(1);
maxx = xlims(2);
miny = ylims(1);
maxy = ylims(2);
origDist.bounds = [minx maxx miny maxy];

% FIRST PASS
close(tmpfig)
tmpfig1 = figure('Position',[1169  360  701  563]);
ntap = -1;
ptap = 50;
h(1,1:3) = [0    ntap 0];
h(2,1:3) = [ntap ptap ntap];
h(3,1:3) = [0    ntap 0];
Y = filter2( h, plotdat );


% % % not an edge candidate
% % idx = find( Y < -300 | Y > -255);
% % Y(idx) = 0;
% % % edge candidate
% % idx = find( Y >= -300 & Y <= -255 );
% % Y(idx) = 1;

% not an edge line
idx = find( Y ~= -255 );
Y(idx) = 0;
% is an edge line
idx = find( Y == -255 );
Y(idx) = 1;
imagesc(Y)
colormap([0 0 0; 1 1 1])


% GET AXES BOUNDARY (PIXELS)
tmpfig2 = figure('position',[127  169  994  754])
h = [0 -1 0; -1 2 -1; 0 -1 0];
Ytemp = filter2( h, plotdat );
Ytemp = Ytemp(2:end-1,2:end-1);
% % not an edge line
% idx = find( Ytemp ~= -255 );
% Ytemp(idx) = 0;
% % is an edge line
% idx = find( Ytemp == -255 );
% Ytemp(idx) = 1;

% not an edge candidate
idx = find( Ytemp < -300 | Ytemp > -255);
Ytemp(idx) = 0;
% edge candidate
idx = find( Ytemp >= -300 & Ytemp <= -255 );
Ytemp(idx) = 1;
imagesc(Ytemp)
colormap([0 0 0; 1 1 1])

[row,col] = find(Ytemp,1,'first');
point1 = [col+2 row];
[row,col] = find(Ytemp,1,'last');
point2 = [col+2 row];

minx = round( min(point1(1,1),point2(1,1)) );
maxx = round( max(point1(1,1),point2(1,1)) );
miny = round( min(point1(1,2),point2(1,2)) );
maxy = round( max(point1(1,2),point2(1,2)) );

pixelBounds = [minx maxx miny maxy];
% reduce plotdat to just the plot area
Y = Y([miny:maxy],[minx:maxx]);

close(tmpfig1); close(tmpfig2);
tmpfig = figure('position',[127  169  994  754]);
imagesc(Y);
colormap([0 0 0; 1 1 1]);

% get just the activated points from Y (the edge points)
[idx1,idx2] = find(Y == 1);
% boundary = Y(idx1,idx2);

nxPixels = size(Y,1);
nyPixels = size(Y,2);
minx = origDist.bounds(1);
maxx = origDist.bounds(2);
miny = origDist.bounds(3);
maxy = origDist.bounds(4);

origDist.drStep = (maxy - miny) ./ nyPixels;
origDist.crStep = (maxx - minx) ./ nxPixels;

dr = idx1 * origDist.drStep + miny;
cr = idx2 * origDist.crStep + minx;

close(tmpfig)
end







