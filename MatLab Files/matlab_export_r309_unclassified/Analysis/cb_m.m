function cb_m(varargin)

    linecolor = 'b';
    tb_facealpha = 0.7;
    
    pick_pts = ginput(2);
    x1 = pick_pts(1,1);
    x2 = pick_pts(2,1);
    y1 = pick_pts(1,2);
    y2 = pick_pts(2,2);
    
    m = (y2 - y1) / (x2 - x1);
    
%     idx_data = find( (1e6*t_on >= x1) & (1e6*t_on <= x2) );
%     y_data = angle(refsig_on(idx_data));
%     x_data = t_on(idx_data);
%     idx_n0 = find(y_data ~= 0);
%     y_data = y_data(idx_n0);
%     x_data = x_data(idx_n0);
%     pres = polyfit(x_data,y_data,1);
%     m = pres(1);
%     b = pres(2);
    line_x = [x1 x2];
    line_y = [y1 y2];
    hLine = line(line_x,line_y);
    set(hLine,'Color',linecolor,'LineWidth',3);
    
    % get X-axis units
    xlab = get(gca,'XLabel');
    if(ishandle(xlab))
        xlab = get(xlab,'String');
    end
    idx1 = strfind(xlab,'('); idx2 = strfind(xlab,')');
    if(isempty(idx1) | isempty(idx2))
        idx1 = strfind(xlab,'['); idx2 = strfind(xlab,']');
    end
    if(isempty(idx1) | isempty(idx2))
        units1 = xlab;
    else
        units1 = xlab(idx1+1:idx2-1);
    end
    
    % get Y-axis units
    ylab = get(gca,'YLabel');
    if(ishandle(ylab))
        ylab = get(ylab,'String');
    end
    idx1 = strfind(ylab,'('); idx2 = strfind(ylab,')');
    if(isempty(idx1) | isempty(idx2))
        idx1 = strfind(ylab,'['); idx2 = strfind(ylab,']');
    end
    if(isempty(idx1) | isempty(idx2))
        units2 = ylab;
    else
        units2 = ylab(idx1+1:idx2-1);
    end
    
    hold on;
    ud.hTextBG = patch([0 1 1 0],[0 0 1 1],'w','Visible','off','FaceAlpha',tb_facealpha,'EdgeColor','none');
    hold off;
    hText = text(0.5*(x1+x2),0.5*(y1+y2),[sprintf(' %4.4g ',m) units2 ' / ' units1]);
    set(hText,'Color',linecolor,'FontWeight','bold','BackgroundColor','none',...
        'HorizontalAlignment','left')
    tbpos = get(hText,'Extent');
    tbBGposx = [tbpos(1) tbpos(1)+tbpos(3) tbpos(1)+tbpos(3) tbpos(1)];
    tbBGposy = [tbpos(2) tbpos(2)          tbpos(2)+tbpos(4) tbpos(2)+tbpos(4)];
    set(ud.hTextBG,'XData',tbBGposx,'YData',tbBGposy,'Visible','on')
    
    
    
    
    
    % link up the line and textbox
    if(isprop(gcf,'WindowButtonMotionFcn'))
        ud.textbox_offset = [1 1];
        ud.hLine = hLine;
        ud.hText = hText;
        ud.units1 = units1;
        ud.units2 = units2;
        ud.hax = gca;
        ud.selected = 0;
        ud.hp_endpts = [0 0];
        hold on;
        endpt_size = 10;
        ud.hp_endpts(1) = plot(line_x(1),line_y(1),'o', ...
            'Color',linecolor,'MarkerFaceColor','none',...
            'Visible','off','MarkerSize',endpt_size,'LineWidth', 2, ...
            'ButtonDownFcn',@endpt_click);
        ud.hp_endpts(2) = plot(line_x(2),line_y(2),'o', ...
            'Color',linecolor,'MarkerFaceColor','none',...
            'Visible','off','MarkerSize',endpt_size, 'LineWidth', 2, ...
            'ButtonDownFcn',@endpt_click);
        ud.hp_ctrpt = plot(0.5*sum(line_x),0.5*sum(line_y),'sq', ...
            'Color',linecolor,'MarkerFaceColor','none', ...
            'Visible','off','MarkerSize',endpt_size, 'LineWidth', 2, ...
            'ButtonDownFcn',@ctrpt_click);
        hold off;
        
        % compute the offset of the textbox
        compute_textbox_offset_fig;

        % create a ButtonDown callback on the line object
        set(hLine,'ButtonDownFcn',@line_click)
        set(hText,'ButtonDownFcn',@textbox_click)

        % set destructors on everything created here
        set(ud.hLine,'DeleteFcn',@cleanup)
        set(ud.hText,'DeleteFcn',@cleanup)
        set(ud.hp_endpts,'DeleteFcn',@cleanup)
        set(ud.hp_ctrpt,'DeleteFcn',@cleanup)
    end
    
    % =====================================================================
    % =====================================================================
    % CALLBACKS AND MISC FUNCTIONS
    % =====================================================================
    % =====================================================================
    function m = compute_slope()
        [xarr,yarr] = get_current_pos;
        x1 = xarr(1);
        x2 = xarr(2);
        y1 = yarr(1);
        y2 = yarr(2);
        m = (y2 - y1) / (x2 - x1);
    end

    function endpt_move(varargin)
        hp = gco;
        hax = get(hp,'Parent');
        cpt = get(hax,'CurrentPoint');
        set(hp,'xdata',cpt(1,1),'ydata',cpt(1,2))
        update_ctrpt
        update_line
    end

    function update_textbox(varargin)
        [xarr,yarr] = get_current_pos;
        [cx,cy]     = get_current_center_pos;
        m = compute_slope;
        c_fig = pointInAxesToFigurePosition([cx cy],gca);
        cx_fig = c_fig(1);
        cy_fig = c_fig(2);
        newc_fig = [ud.textbox_offset_fig(1)+cx_fig ud.textbox_offset_fig(2)+cy_fig];
        newc_axe = pointInFigureToAxesPosition(newc_fig,gca);
        
        set(ud.hText,'Position',[newc_axe(1) newc_axe(2)], ...
            'String',[sprintf(' %4.4g ',m) units2 ' / ' units1]);
        update_textbox_background;
    end

    %======================================================================
    % BELOW ARE GENERIC FUNCTIONS THAT SHOULD BE COMMON FOR cb_dx.m,
    % cb_dy.m, cb_m.m, cb_frq.m
    %======================================================================
    function line_click(varargin)
        if(ud.selected == 0)
            ud.selected = 1;
            set(ud.hp_endpts,'Visible','on');
            set(ud.hp_ctrpt,'Visible','on');
        else
            ud.selected = 0;
            set(ud.hp_endpts,'Visible','off');
            set(ud.hp_ctrpt,'Visible','off');
        end
    end

    function endpt_click(varargin)
        h_endpt = gco;
        hax  = get(h_endpt,'Parent');
        hfig = get(hax,'Parent');
        % set callbacks for WindowButton callbacks interacting with this
        % point
        set(hfig,'WindowButtonMotionFcn',@endpt_move)
        set(hfig,'WindowButtonUpFcn',@button_release)
    end

    function ctrpt_click(varargin)
        h_endpt = gco;
        hax  = get(h_endpt,'Parent');
        hfig = get(hax,'Parent');
        % set callbacks for WindowButton callbacks interacting with this
        % point
        set(hfig,'WindowButtonMotionFcn',@ctrpt_move)
        set(hfig,'WindowButtonUpFcn',@button_release)
    end

    function textbox_click(varargin)
        h_tb = gco;
        hax  = get(h_tb,'Parent');
        hfig = get(hax,'Parent');
        % set callbacks for WindowButton callbacks interacting with this
        % point
        set(hfig,'WindowButtonMotionFcn',@textbox_move)
        set(hfig,'WindowButtonUpFcn',@button_release)
    end

    function button_release(varargin)
        set(gcf,'WindowButtonMotionFcn',[])
        set(gcf,'WindowButtonUpFcn',[])
    end

    function ctrpt_move(varargin)
        hp = gco;
        hax = get(hp,'Parent');
        cpt = get(hax,'CurrentPoint');
        set(hp,'xdata',cpt(1,1),'ydata',cpt(1,2))
        update_location
    end

    function textbox_move(varargin)
        htb = gco;
        hax = get(htb,'Parent');
        cpt = get(hax,'CurrentPoint');
        
        set(htb,'Position',[cpt(1,1) cpt(1,2)])
        
        compute_textbox_offset_fig;
        
        update_textbox_background;
    end

    function update_line()
        [xarr,yarr] = get_current_pos;
        set(ud.hLine,'XData',xarr,'YData',yarr);
        update_textbox;
    end

    function update_location()
        
        [cx,cy]     = get_current_center_pos();
        [xarr,yarr] = get_current_pos;
        
        dx = diff(xarr);
        dy = diff(yarr);
        
        new_xarr = cx + 0.5*dx*[-1 1];
        new_yarr = cy + 0.5*dy*[-1 1];
        
        set(ud.hp_endpts(1),'XData',new_xarr(1),'YData',new_yarr(1))
        set(ud.hp_endpts(2),'XData',new_xarr(2),'YData',new_yarr(2))
        update_line;
        update_textbox;
    end

    function update_textbox_background(varargin)
        tbpos = get(ud.hText,'Extent');
        % tbpos: x, y, wx, wy
        %   y happens to be the corner
        if(strcmpi(get(ud.hax,'XDir'),'reverse'))
            tbBGposx = [tbpos(1) tbpos(1)-tbpos(3) tbpos(1)-tbpos(3) tbpos(1)];
        else
            tbBGposx = [tbpos(1) tbpos(1)+tbpos(3) tbpos(1)+tbpos(3) tbpos(1)];
        end
        if(strcmpi(get(ud.hax,'YDir'),'reverse'))
            tbBGposy = [tbpos(2) tbpos(2)          tbpos(2)-tbpos(4) tbpos(2)-tbpos(4)];
        else
            tbBGposy = [tbpos(2) tbpos(2)          tbpos(2)+tbpos(4) tbpos(2)+tbpos(4)];
        end
        
        set(ud.hTextBG,'XData',tbBGposx,'YData',tbBGposy,'Visible','on')
    end

    function update_ctrpt()
        [xarr,yarr] = get_current_pos;
        cx = 0.5*sum(xarr);
        cy = 0.5*sum(yarr);
        set(ud.hp_ctrpt,'xdata',cx,'ydata',cy)
    end



    function [xarr,yarr] = get_current_pos()
        % get current point positions
        x1 = get(ud.hp_endpts(1),'XData');
        x2 = get(ud.hp_endpts(2),'XData');
        y1 = get(ud.hp_endpts(1),'YData');
        y2 = get(ud.hp_endpts(2),'YData');
        xarr = [x1 x2];
        yarr = [y1 y2];
    end

    function [cx,cy] = get_current_center_pos()
        cx = get(ud.hp_ctrpt,'XData');
        cy = get(ud.hp_ctrpt,'YData');
    end

    function compute_textbox_offset_fig()
        [cx,cy] = get_current_center_pos;
        c_fig  = pointInAxesToFigurePosition([cx cy],gca);
        cx_fig = c_fig(1);
        cy_fig = c_fig(2);
        tbpos = get(ud.hText,'Position');
        t_fig = pointInAxesToFigurePosition(tbpos(1:2), gca);
        tx_fig = t_fig(1);
        ty_fig = t_fig(2);
        ud.textbox_offset_fig = [tx_fig-cx_fig  ty_fig-cy_fig];
    end

    function cleanup(varargin)
        set(ud.hLine,'DeleteFcn',[])
        set(ud.hText,'DeleteFcn',[])
        set(ud.hTextBG,'DeleteFcn',[])
        set(ud.hp_endpts,'DeleteFcn',[])
        set(ud.hp_ctrpt,'DeleteFcn',[])
        delete(ud.hLine)
        delete(ud.hText)
        delete(ud.hTextBG)
        delete(ud.hp_endpts)
        delete(ud.hp_ctrpt)
    end


end