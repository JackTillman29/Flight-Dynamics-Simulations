function cb_grab_slice_vertical(varargin)
    pick_pts = ginput(1);
    x1 = pick_pts(1,1);
    y1 = pick_pts(1,2);
    
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
    
    % get title string from axes
    titlestr  = get(get(haxorig,'Title'),'String');
    xlabelstr = get(get(haxorig,'XLabel'),'String');
    ylabelstr = get(get(haxorig,'YLabel'),'String');
    try
        hcb = findobj(get(gcf,'Children'),'Type','ColorBar');
        clabelstr = get(get(hcb,'YLabel'),'String');
    catch
        clabelstr = '';
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
    
    
    if(strcmpi(get(datobj,'Type'),'line'))
        % if object is a line object (when expecting a 2D data object),
        % then go to the axes object, and look for a 2D data object child
        %  THIS IS A COMMON ERROR WHERE A LINE DRAWN ON THE PLOT FROM
        %  PREVIOUS CUTS MADE ON THE IMAGE WOULD BE DETECTED.
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
    [xdat,ydat,zq] = get_data_vert(datobj,x1);
    
    xline = [x1 x1];
    yline = [ylims(1) ylims(2)];
    hLine = line(xline,yline);
    set(hLine,'Color','r','LineWidth',3);
    
    hfig_slice = figure;
    haxnew = axes;
    add_print_callbacks;
    add_analysis_callbacks;
    hp_slice_data = plot(haxnew,ydat,zq);
   
    % because this is a vertical slice, put original plot's YLABEL on this
    % cut plot's XLABEL
    xlabel(haxnew,ylabelstr);
    ylabel(haxnew,unit_label);
    title(haxnew,titlestr);
    
    % apply original plot's yticks to this plot's XTICKS
    set(haxnew,'xtick',yticks);
    % apply original plot's ylim to this plot's XLIM
    set(haxnew,'xlim',ylims);
    % apply original plot's caxis to this plot's YLIM
    set(haxnew,'ylim',clims);
    
    try
        set(hLine,'ButtonDownFcn',@button_down_line);
        set(hfig_datobj,'WindowButtonUpFcn',@button_up_line_vert);
        set(hfig_datobj,'WindowButtonMotionFcn',@move_line_vert);
    catch
        
    end

    % link the drawn figure to this line
    ud_slice_line.hfig_slice = hfig_slice;
    ud_slice_line.hp_slice_data = hp_slice_data;
    ud_slice_line.hfig_datobj = hfig_datobj;
    set(hLine,'UserData',ud_slice_line);
    set(hLine,'Tag','slice_line_vert')

    function [xdat,ydat,zq] = get_data_vert(datobj,x1)
        
        switch get(datobj,'Type')
            case 'image'
                xdat  = get(datobj,'XData');
                ydat  = get(datobj,'YData');
                z     = get(datobj,'CData');
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
                if( (ix1 < 1) )
                    zq = z(:,1);
                    xdat = min(xdat);
                elseif( ix1 > size(z,2) )
                    zq = z(:,end);
                    xdat = max(xdat);
                else
                    zq = z(:,ix1);
                end

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
                if( (ix1 < 1) )
                    zq = z(1,:);
                    xdat = min(xdat);
                elseif( ix1 > size(z,1) )
                    zq = z(end,:);
                    xdat = max(xdat);
                else
                    zq = z(ix1,:);
                end
        end
    end

    function button_down_line(varargin)
        if(strcmpi(get(gco,'Type'),'line'))
            ud = get(gco,'UserData');
            % if the line that was clicked is the slice line
            if(isfield(ud,'hfig_slice'))
                ud.slice_clicked = 1;
                set(gco,'UserData',ud)
                set(hfig_datobj,'WindowButtonUpFcn',@button_up_line_vert)
                set(hfig_datobj,'WindowButtonMotionFcn',@move_line_vert)
            end
        end
    end

    function button_up_line_vert(varargin)
        % find object in the current figure with tag "slice_line_vert"
        hslice = findobj(get(gcf,'children'),'Tag','slice_line_vert');
        for k = 1:length(hslice)
            ud = get(hslice(k),'UserData');
            % if the line that was clicked is the slice line
            if(isfield(ud,'slice_clicked'))
                % if multiple slices have been clicked
                if(ud.slice_clicked == 1)
                    ud.slice_clicked = 0;
                    set(hslice(k),'UserData',ud)
                end
            end
        end
    end

    function move_line_vert(varargin)
        hfig = varargin{1};
        hslice = findobj(get(gcf,'children'),'Tag','slice_line_vert');
        for k = 1:length(hslice)
            ud = get(hslice(k),'UserData');
            if(isfield(ud,'slice_clicked'))
                if(ud.slice_clicked == 1)
                    hax = get(hslice(k),'Parent');
                    cpt = get(hax,'CurrentPoint');
                    new_xdat = cpt(1,1);
                    
                    [xdat,ydat,zq] = get_data_vert(datobj,new_xdat);
                    set(ud.hp_slice_data,'xdata',ydat,'ydata',zq)
                    if(length(xdat) == 1)
                        set(hslice(k),'xdata',xdat*ones(1,2));
                    else
                        set(hslice(k),'xdata',new_xdat*ones(1,2));
                    end
%                     disp('here!')
%                     varargin{1}
%                     varargin{2}
%                     length(varargin)
                end
            end
        end
    end





end
