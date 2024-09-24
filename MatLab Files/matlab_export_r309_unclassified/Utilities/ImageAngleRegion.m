function [maxy,avgy,max_xloc,max_yloc,max_rng,meanval,stdval] = ImageAngleRegion(xgrid,ygrid,cdata,adata,rng,origin,vertang,angwidth,plotflag)
% Purpose of this function is to compute an average aspect value based upon
% the inputs provided
% K. Sawmiller, BAH, 2014
%
% Things we need: 
% 1 - Image CData,xgrid,ygird
% 2 - Initial range of interest
% 3 - 0,0 point
% 4 - Angle from vertical
% 5 - Aspect width from (4) in degrees
%
% Returns [maxy,avgy,max_xloc,max_yloc,max_rng,meanval,stdval]
% maxy is maximum value in the region of interest
% avgy is average value in the region of interest
% max_xloc is maximum value x position
% max_yloc is maximum value y position
% max_rng is maximum value range from origin
% meanval is average value along rng input arc
% stdval is standard deviation along rng input arc
if(~exist('plotflag'))
    plotflag = 0;
end

%adata_in = adata;
%adata = 0.0 * adata_in + 1;



[XGRID,YGRID] = meshgrid(xgrid-origin(1),ygrid-origin(2));
vert_angle_grid = atan2(XGRID,YGRID);

RGRID = sqrt(XGRID.^2 + YGRID.^2);

ang_test_min = vertang-angwidth/2.0;
ang_test_max = vertang+angwidth/2.0;

% argument check. if both angles are > 180, adjust to -angles
if(ang_test_min > pi)
    ang_test_min = ang_test_min - 2*pi;
    ang_test_max = ang_test_max - 2*pi;
end

cdata2 = ( ...
    adata == 1 & ...                                % within detection range
    RGRID >= rng & ...                              % outside of range of interest
    (vert_angle_grid >= (ang_test_min) ) & ... % within angle cut of interest
    vert_angle_grid <= (ang_test_max) );

% odd case when +/- 180 is included in data set
if(ang_test_max >= pi & ang_test_min <= pi)
            cdata2 = cdata2 + ...
            ( adata == 1 & ...                                % within detection range
              RGRID >= rng & ...                              % outside of range of interest
              (vert_angle_grid <= (ang_test_max-2*pi) ));

end



xi = find(cdata2 ~=0);
azi = vert_angle_grid(xi);

[maxy,maxi]=max(cdata(xi));

temp=XGRID(xi);
max_xloc = temp(maxi);
temp=YGRID(xi);
max_yloc = temp(maxi);
temp=RGRID(xi);
max_rng  = temp(maxi);

% average calculation
avgy = mean(cdata(xi));

% calulations along minimum range arc
cdata3 = cdata2 & RGRID < (rng + 2*mean(diff(xgrid)));

rng_arc_x = 1:size(cdata3,2);
rng_arc_y = zeros(1,size(cdata3,2));


for k = 1 : length(rng_arc_x)
    a=find(cdata3(:,k) == 1, 1, 'first');
    if(isempty(a))
        rng_arc_y(k) = nan;
        rng_arc_x(k) = nan;
    else
        rng_arc_y(k) = a;
    end
end

% Create index i,j arrays for curve
rng_arc_x = rng_arc_x(find(~isnan(rng_arc_x)));
rng_arc_y = rng_arc_y(find(~isnan(rng_arc_y)));

values = zeros(1,length(rng_arc_x));
for k = 1 : length(rng_arc_x)
    values(k) = cdata(rng_arc_y(k),rng_arc_x(k));
end

meanval = mean(values);
stdval  = std(values);


if(isempty(maxy))
    maxy = nan;
end
if(isempty(max_rng))
    max_rng = nan;
end



for k = 1 : length(rng_arc_y)
    if(~isnan(rng_arc_y(k)))
        rng_arc_y(k) = ygrid(rng_arc_y(k));
    else
        rng_arc_y(k) = nan;
    end
end
rng_arc_x = xgrid(rng_arc_x);


% debug region
if(nargout == 0 || plotflag == 1)
    figure;
    set(gcf,'Color','k');
    colormap(jet(1024));
    hi = imagesc(xgrid,ygrid,cdata);
    set(hi,'AlphaData',0.5*(0.75*adata+1.25*cdata2));
    set(gca,'YDir','normal','Color','k','XColor',0.4*[1 1 1],'YColor',0.4*[1 1 1]);
    hold on;
    plot(origin(1),origin(2),'wx');
    axis equal;
    grid on;
    xlabel('Crossrange');ylabel('Downrange');hc=colorbar;
    set(hc,'YColor',0.4*[1 1 1],'XColor',0.4*[1 1 1]);
    hold on;
    plot(origin(1)+max_xloc,origin(2)+max_yloc,'wo');
    %ht=text(origin(1)+max_xloc,origin(2)+max_yloc,['  R=' num2str(max_rng,'%5.1f') ', MAX=' num2str(maxy,'%5.1f')]);
    title(['Region Average=  ' num2str(avgy)],'Color',0.4*[1 1 1]);
%    set(ht,'Color',0.8*[1 1 1],'FontWeight','bold');
    
    set(hi,'ButtonDownFcn','upd_ImageAngleRegion');
    set(gcf,'DoubleBuffer','on');
    
    % min range arc
    line(rng_arc_x,rng_arc_y,'Color','w','LineWidth',2);
    try
        ht=text(rng_arc_x(end),rng_arc_y(end), ...
            {['  AVG=' num2str(meanval,'%5.1f') ', STD=' num2str(stdval,'%5.1f')], ...
            ['  R=' num2str(max_rng,'%5.1f') ', MAX=' num2str(maxy,'%5.1f')]});
        
        set(ht,'Color',0.8*[1 1 1],'FontWeight','bold');
    catch
        disp('Could not create text marker.');
    end
    
end

end