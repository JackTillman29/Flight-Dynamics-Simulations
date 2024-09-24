function cb_grab_slice_slant(varargin)
    pick_pts = ginput(2);
    x1 = pick_pts(1,1);
    x2 = pick_pts(2,1);
    y1 = pick_pts(1,2);
    y2 = pick_pts(2,2);
    
    interpolation_method = 'nearest';
%     interpolation_method = 'linear';
%     interpolation_method = 'cubic';
    
    haxorig = gca;
    hchild = get(haxorig,'children');
    
    % find the first image or surface child of this axes
    %     children handles are ordered in a "pile-on" format, where the
    %     oldest child (first object created in axes) is on the bottom, and
    %     the youngest child (last object created in axes) is on the top,
    %     so flipud(hchild) will order it the opposite.
    hchild = flipud(hchild);
    for k = 1:length(hchild)
        switch get(hchild(k),'Type')
            case 'image'
                datobj = hchild(k);
            case 'surface'
                datobj = hchild(k);
        end
    end
    
    % get data in the axes
    xlims = get(haxorig,'XLim');
    ylims = get(haxorig,'YLim');
    
    % compute slope
    m = (y2 - y1) / (x2 - x1);
    %[x1 x2 y1 y2]
    
    % get title string from axes
    titlestr  = get(get(haxorig,'Title'),'String');
    xlabelstr = get(get(haxorig,'XLabel'),'String');
    ylabelstr = get(get(haxorig,'YLabel'),'String');
    try
        hcb = findobj(get(gcf,'Children'),'Type','ColorBar');
        clabelstr = get(get(hcb,'YLabel'),'String');
    catch
        try
%             hcb = 
        catch
            clabelstr = '';
        end
    end
    clims = caxis(haxorig);
    xticks = get(haxorig,'XTick');
    yticks = get(haxorig,'YTick');
    if(isempty(titlestr))
        titlestr = 'Image Data';
    end
    if(isempty(xlabelstr))
        xlabelstr = 'X';
    end
    if(isempty(ylabelstr))
        ylabelstr = 'Y';
    end
    
    get(datobj,'Type')
    if(strcmpi(get(datobj,'Type'),'line'))
        % if object is a line object (when expecting a 2D data object),
        % then go to the axes object, and look for a 2D data object child
        hchild = get(haxorig,'Children');
        datobj = findobj(hchild,'Type','image');
        if(isempty(datobj))
            datobj = findobj(hchild,'Type','surface');
            if(isempty(datobj))
                error('couldn''t find a valid data object on the figure')
            end
            % if there are multiple image objects on the current axes
            % (how/why?), then get the FIRST one.
            if(length(datobj) > 1)
                datobj = datobj(1);
            end
        else
            % if there are multiple image objects on the current axes
            % (how/why?), then get the FIRST one.
            if(length(datobj) > 1)
                datobj = datobj(1);
            end
        end
    end
    
    % see if axes containing the image/surface has a colorbar, then get
    % any label attached to it
    hax_datobj = get(datobj,'Parent');
    hfig_datobj = get(hax_datobj,'Parent');
    try
        hcb_datobj = get(hax_datobj,'Colorbar');
        unit_label = get(get(hcb_datobj,'Label'),'String');
    catch
        hcb_datobj = findobj(hfig_datobj,'tag','colorbar');
        unit_label = get(get(hcb_datobj,'ylabel'),'string');
    end
    
    
    % capture the indices into the data for these start and end points
    switch get(datobj,'Type')
        case 'image'
            xdat = get(datobj,'XData');
            ydat = get(datobj,'YData');
            z = get(datobj,'CData');
            zsize = size(z);
            if(length(xdat) ~= size(z,2))
                % image x data is just the LIMITS of the axes
                xdat = linspace(xdat(1),xdat(2),size(z,2));
            end
            if(length(ydat) ~= size(z,1))
                % image x data is just the LIMITS of the axes
                ydat = linspace(ydat(1),ydat(2),size(z,1));
            end
            nx = zsize(2);
            ny = zsize(1);
            % convert x1 and x2 into INDICES
            minmax_xdat = [min(xdat) max(xdat)];
            minmax_ydat = [min(ydat) max(ydat)];
            dxdat = abs(diff(minmax_xdat));
            dydat = abs(diff(minmax_ydat));
            ix1 = ((x1-min(xdat)) / dxdat)*nx;
            ix2 = ((x2-min(xdat)) / dxdat)*nx;
            iy1 = ((y1-min(ydat)) / dydat)*ny;
            iy2 = ((y2-min(ydat)) / dydat)*ny;
            ixstep = sign(ix2-ix1);
            iystep = sign(iy2-iy1);
            ix = ix1:ixstep:ix2;
            iy = iy1:iystep:iy2;
            mx = length(ix) % # of points along x side
            my = length(iy) % # of points along y side
            mz = round(sqrt(mx^2 + my^2)) % # of points along hypotenuse
            
            ixq = round(linspace(ix1,ix2,mz));
            iyq = round(linspace(iy1,iy2,mz));
            
            [IX,IY] = meshgrid(1:nx,1:ny);
            zq = interp2(IX,IY,z,ixq,iyq,interpolation_method);
            
        case 'surface'
            % todo
            xdat = get(datobj,'XData');
            ydat = get(datobj,'YData');
            z = get(datobj,'ZData');
            zsize = size(z);
            
            if(length(ydat) ~= size(z,1))
                % image x data is just the LIMITS of the axes
                ydat = linspace(ydat(1),ydat(2),size(z,1));
            end
            
            nx = zsize(2);
            ny = zsize(1);
            % convert x1 and x2 into INDICES
            minmax_xdat = [min(xdat) max(xdat)];
            minmax_ydat = [min(ydat) max(ydat)];
            dxdat = abs(diff(minmax_xdat));
            dydat = abs(diff(minmax_ydat));
            ix1 = round(((x1-min(xdat)) / dxdat)*nx);
            zq = z(:,ix1);
    end
    
%     if(~isempty(xlab.String))
%         xidx1 = strfind(xlab,'(');
%         xidx2 = strfind(xlab,')');
%         xlab = xlab(xidx1:xidx2);
%     end
%     if(~isempty(xlab.String))
%         yidx1 = strfind(ylab,'(');
%         yidx2 = strfind(ylab,')');
%         ylab = ylab(yidx1:yidx2);
%     end
    
    
%     if(strcmpi(xlab,'s'))
%         timeunits = 1;
%     elseif(strcmpi(xlab,'ms'))
%         timeunits = 1e-3;
%     elseif(strcmpi(xlab,'\mus'))
%         timeunits = 1e-6;
%     elseif(strcmpi(xlab,'ns'))
%         timeunits = 1e-9;
%     end
%     timestr = xlab;
%     
%     if(strcmpi(ylab,'Hz'))
%         frequnits = 1;
%     elseif(strcmpi(ylab,'kHz'))
%         frequnits = 1e3;
%     elseif(strcmpi(ylab,'MHz'))
%         frequnits = 1e6;
%     elseif(strcmpi(ylab,'GHz'))
%         frequnits = 1e9;
%     end
%     freqstr = ylab;
    
    xline = [x1 x2];
    yline = [y1 y2];
    hLine = line(xline,yline);
    set(hLine,'Color','r','LineWidth',3);
    
    
    hfignew = figure;
    haxnew = axes;
    add_print_callbacks;
    add_analysis_callbacks;
    plot(haxnew,zq);
    % because this is a slanted slice, put a label for "indices" on XLABEL
    xlabel('Indices');
    ylabel(unit_label);
    title(titlestr);
    dat.x = linspace(x1,x2,length(zq));
    dat.y = linspace(y1,y2,length(zq));
    dat.z = zq;
    dat.xlab = xlabelstr;
    dat.ylab = ylabelstr;
    dat.zlab = clabelstr;
    set(hfignew,'UserData',dat);
    ylim(haxnew,clims);
    
    dcm_obj = datacursormode(hfignew);
    set(dcm_obj,'UpdateFcn',@dcm_update_fcn);
    
    
%     xticks  = get(hax,'xtick')
%     xticks(find(xticks == 0)) = 1;
%     xticks(find(xticks > length(zq))) = length(zq);
%     
%     newxticks = cell(1,length(xticks));
%     for k = 1:length(xticks)
%         newxticks{k} = [num2str([dat.x(xticks(k)); dat.y(xticks(k))]) [char(10);' ']];
%     end
%     set(hax,'xtick',xticks,'XTickLabel',newxticks)
    
    
%     hText = text(x2,y2,sprintf([' BW = %4.1f ',freqstr],1/m))
%     set(hText,'Color','b','FontWeight','bold','BackgroundColor','w')
    
    
%     t_on = evalin('base','t_on');
%     refsig_on = evalin('base','refsig_on');
%     
%     idx_data = find( (1e6*t_on >= x1) & (1e6*t_on <= x2) );
%     y_data = angle(refsig_on(idx_data));
%     x_data = t_on(idx_data);
%     pres = polyfit(x_data,y_data,1);
%     m = pres(1);
%     b = pres(2);
%     line_x = [x1 x2];
%     line_y = [ ...
%         1e-6*m*x1 + b ...
%         1e-6*m*x2 + b];
%     hLine = line(line_x,line_y);
%     set(hLine,'Color','b','LineWidth',3);
%     hText = text(x2,double(line_y(end)),sprintf(' %4.1f Hz',m));
%     set(hText,'Color','b','FontWeight','bold','BackgroundColor','w')

function txt = dcm_update_fcn(varargin)

var1 = varargin{1};
event_obj = varargin{2};


hax = get(var1,'Parent');
hfig = get(hax,'Parent');

% get original xy array data from figure's "UserData" field
dat = get(hfig,'UserData');

pos = get(event_obj,'Position');
idx = pos(1);

txt = {'Data from image at:'; ...
    [dat.xlab ': ' num2str(dat.x(idx))]; ...
    [dat.ylab ': ' num2str(dat.y(idx))]; ...
    [dat.zlab ': ' num2str(dat.z(idx))]};

end

end
