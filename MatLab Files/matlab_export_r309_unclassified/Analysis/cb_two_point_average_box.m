function cb_two_point_average_box(varargin)
    
    patchcolor = 'w';
    patchfacecolor = 'none';
    tb_facealpha = 0.7;
    
    pt_facecolor = 'none';
    pt_edgecolor = 'w';
    
    x=ginput(2);
    
    ud.hax = gca;
    ud.hfig = get(ud.hax,'Parent');
    
    ud.hi=findobj('Type','image','Parent',gca);
    
    minx = min(x(:,1));
    maxx = max(x(:,1));
    miny = min(x(:,2));
    maxy = max(x(:,2));
    midx = 0.5*sum(x(:,1));
    midy = 0.5*sum(x(:,2));
    edgeL = [minx midy];
    edgeR = [maxx midy];
    edgeB = [midx miny];
    edgeT = [midx maxy];
    
    cornerUL = [minx maxy];
    cornerUR = [maxx maxy];
    cornerBL = [minx miny];
    cornerBR = [maxx miny];
    
    %=============
    fetch_data([minx maxx],[miny maxy])
    %=============

    % if linear
    %avg1 = mean(c_window(:));

    % if log
    %avg2 = 10*log10(mean(10.^(0.1*c_window(:))));
    %disp('lin,log');
    %disp([avg1,avg2])
    %msgbox(sprintf('Linear Average: %18.6e\nLog Average: %8.3f',avg1,avg2),'Average');
    ud.hpatch = patch( ...
        [x(1,1) x(2,1) x(2,1) x(1,1)], ...
        [x(1,2) x(1,2) x(2,2) x(2,2)], ...
        'r');
    set(ud.hpatch,'EdgeColor',patchcolor,'FaceColor',patchfacecolor);
    ud.data_kind = questdlg('Linear or Log (dB) Data?','Data Question','Linear','Log','Log');
    
    %=============
    compute_stats;
    %=============
    
    midx = mean(x(:,1));
    midy = mean(x(:,2));
    
    ud.htextBG = patch([0 1 1 0],[0 0 1 1],'w','Visible','off','FaceAlpha',tb_facealpha,'EdgeColor','none');
    ud.htext = text(midx,midy,ud.stat_str);
    set(ud.htext,'Color','k','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
%     tbpos = get(ud.htext,'Extent')
%     tbBGposx = [tbpos(1) tbpos(1)+tbpos(3) tbpos(1)+tbpos(3) tbpos(1)]
%     tbBGposy = [tbpos(2) tbpos(2)          tbpos(2)+tbpos(4) tbpos(2)+tbpos(4)]
%     set(ud.htextBG,'XData',tbBGposx,'YData',tbBGposy,'Visible','on')

    hold on;
    ud.hp_edges(1)  = plot(edgeL(1),edgeL(2),'sq');
    ud.hp_edges(2)  = plot(edgeR(1),edgeR(2),'sq');
    ud.hp_edges(3)  = plot(edgeB(1),edgeB(2),'sq');
    ud.hp_edges(4)  = plot(edgeT(1),edgeT(2),'sq');
    ud.hp_ctrpt      = plot(midx,    midy,    'sq');
    
    % corner order:   UL, UR, BL, BR
    ud.hp_corners(1) = plot(cornerUL(1),cornerUL(2),'sq');
    ud.hp_corners(2) = plot(cornerUR(1),cornerUR(2),'sq');
    ud.hp_corners(3) = plot(cornerBL(1),cornerBL(2),'sq');
    ud.hp_corners(4) = plot(cornerBR(1),cornerBR(2),'sq');
%     ud.hp_edges   = [ud.hp_edges_L ud.hp_edges_R ud.hp_edges_B ud.hp_edges_T];
%     ud.hp_corners = [ud.hp_corners_UL ud.hp_corners_UR ud.hp_corners_BL ud.hp_corners_BR];
    set(ud.hp_edges,  'Visible','off')
    set(ud.hp_ctrpt, 'Visible','off')
    set(ud.hp_corners,'Visible','off')
    set(ud.hp_edges,  'MarkerFaceColor',pt_facecolor,'Color',pt_edgecolor,'MarkerSize',10,'LineWidth',2);
    set(ud.hp_ctrpt, 'MarkerFaceColor',pt_facecolor,'Color',pt_edgecolor,'MarkerSize',10,'LineWidth',2);
    set(ud.hp_corners,'MarkerFaceColor',pt_facecolor,'Color',pt_edgecolor,'MarkerSize',10,'LineWidth',2);
    hold off;
    
    %=============
    update_textbox_background;
    compute_textbox_offset_fig;
    set(ud.htextBG,'Visible','on')
    %=============
    
    
    % add callbacks to each graphics object
    
    if exist('addlistener') && isprop(ud.hfig,'WindowButtonMotionFcn')
        ud.selected = 0;
        set(ud.hpatch,    'ButtonDownFcn',@click_patch)
        set(ud.hp_edges,  'ButtonDownFcn',@click_edge)
        set(ud.hp_corners,'ButtonDownFcn',@click_corner)
        set(ud.hp_ctrpt,  'ButtonDownFcn',@click_ctr)
        set(ud.htext,     'ButtonDownFcn',@click_textbox)

        set(ud.hpatch,    'DeleteFcn',@cleanup)
        set(ud.htext,     'DeleteFcn',@cleanup)
        set(ud.hp_edges,  'DeleteFcn',@cleanup)
        set(ud.hp_corners,'DeleteFcn',@cleanup)
%         set(ud.hpatch,'DeleteFcn',@cleanup)
    end



    %======================================================================
    %======================================================================
    %  INTERACTIVE CALLBACKS AND FUNCTIONS
    %======================================================================
    %======================================================================
    % Major idea behind managing the edges and corners:
    %  BLUF: whenever changing the points, always update the edges,
    %  corners, and centerpoints in this ORDER:
    %       center
    %
    %       CORNERS MODIFY EDGES, EDGES MODIFY PATCH
    %       MODIFYING CENTER MODIFIES EDGES, EDGES MODIFY CORNERS
    %
    %       FIRST OBJECT -> SECOND OBJECT -> THIRD OBJECT
    %       CORNER -> EDGE -> CENTER -> PATCH
    %       EDGE -> CORNER -> CENTER -> PATCH
    %       CENTER -> EDGE -> CORNER -> PATCH
    function fetch_data(xextent,yextent)
        if(nargin == 0)
            xextent = [];
            yextent = [];
        end
        
        c = get(ud.hi,'CData');
        
        xdata = get(ud.hi,'XData');
        ydata = get(ud.hi,'YData');
        
        % if length(xdata) or size(xdata,1) does not equal size(cdata,1),
        % then create an equally spaced x-vector
        if(length(xdata) ~= size(c,2))
            % most likely, the xdata is the xlimits
            xlims = xdata
            xdata = linspace(xlims(1),xlims(2),size(c,2));
        end
        if(length(ydata) ~= size(c,1))
            % most likely, the xdata is the xlimits
            ylims = ydata
            ydata = linspace(ylims(1),ylims(2),size(c,1));
        end
        
        xa=find(xdata > xextent(1),1,'first');
        xb=find(xdata < xextent(2),1,'last');
        ya=find(ydata > yextent(1),1,'first');
        yb=find(ydata < yextent(2),1,'last');
        xstep = sign(xb - xa);
        ystep = sign(yb - ya);

        ud.c_window = c(ya:ystep:yb,xa:xstep:xb);
    end

    function compute_stats()
        if(strcmpi(ud.data_kind,'Linear'))
            avgv = mean(ud.c_window(:));
            minv = min(ud.c_window(:));
            maxv = max(ud.c_window(:));
            fmt = 'min: %8.3e\nmax: %8.3e\navg: %8.3e';
        else
            avgv = 10*log10(mean(10.^(0.1*ud.c_window(:))));
            minv = 10*log10(min(10.^(0.1*ud.c_window(:))));
            maxv = 10*log10(max(10.^(0.1*ud.c_window(:))));
            fmt = 'min: %8.3f\nmax: %8.3f\navg: %8.3f';
        end
        ud.stat_str = sprintf(fmt,minv,maxv,avgv);
    end

    function update_patch()
        [xarr,yarr] = get_current_corner_pts;
        set(ud.hpatch,'xdata',xarr,'ydata',yarr)
        fetch_data([min(xarr) max(xarr)],[min(yarr) max(yarr)]);
        compute_stats
        update_textbox_stats
    end


    function [xarr,yarr] = get_current_edge_pts()
        % 
        xarr(1) = get(ud.hp_edges(1),'xdata');
        xarr(2) = get(ud.hp_edges(2),'xdata');
        xarr(3) = get(ud.hp_edges(3),'xdata');
        xarr(4) = get(ud.hp_edges(4),'xdata');
        
        yarr(1) = get(ud.hp_edges(1),'ydata');
        yarr(2) = get(ud.hp_edges(2),'ydata');
        yarr(3) = get(ud.hp_edges(3),'ydata');
        yarr(4) = get(ud.hp_edges(4),'ydata');
    end

    function [xarr,yarr] = get_current_corner_pts()
        % 
        for k = 1:4
            xarr(k) = get(ud.hp_corners(k),'xdata');
            yarr(k) = get(ud.hp_corners(k),'ydata');
        end
%         xarr(2) = get(ud.hp_corners(2),'xdata');
%         xarr(3) = get(ud.hp_corners(3),'xdata');
%         xarr(4) = get(ud.hp_corners(4),'xdata');
%         
%         yarr(1) = get(ud.hp_corners(1),'ydata');
%         yarr(2) = get(ud.hp_corners(2),'ydata');
%         yarr(3) = get(ud.hp_corners(3),'ydata');
%         yarr(4) = get(ud.hp_corners(4),'ydata');
    end

    function [xarr,yarr] = organize_corners(xarr,yarr)
        minx = min(xarr);
        maxx = max(xarr);
        miny = min(yarr);
        maxy = max(yarr);
        % UL, UR, BL, BR
        xarr = [minx maxx maxx minx];
        yarr = [miny miny maxy maxy];
%         corners.UL = [minx maxy];
%         corners.UR = [maxx maxy];
%         corners.BL = [minx miny];
%         corners.BR = [maxx miny];
    end

    function [cx,cy] = get_current_ctr_pt()
        cx = get(ud.hp_ctrpt,'xdata');
        cy = get(ud.hp_ctrpt,'ydata');
    end

    function click_patch(varargin)
        if(ud.selected == 0)
            ud.selected = 1;
            set(ud.hp_edges,'Visible','on');
            set(ud.hp_ctrpt,'Visible','on');
            set(ud.hp_corners,'Visible','on');
        else
            ud.selected = 0;
            set(ud.hp_edges,'Visible','off');
            set(ud.hp_ctrpt,'Visible','off');
            set(ud.hp_corners,'Visible','off');
        end
    end

    function click_edge(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_edge)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
        set(ud.hfig,'Interruptible','off')
    end

    function click_ctr(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_ctr)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
        set(ud.hfig,'Interruptible','off')
    end

    function click_corner(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_corner)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
        set(ud.hfig,'Interruptible','off')
    end

    function click_textbox(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_textbox)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
        set(ud.hfig,'Interruptible','off')
    end

    function move_edge(varargin)
        hp = gco;
        cpt = get(ud.hax,'CurrentPoint');
        
        ex = get(hp,'xdata');
        ey = get(hp,'ydata');
        cx = get(ud.hp_ctrpt,'xdata');
        cy = get(ud.hp_ctrpt,'ydata');
        % determine if this edge is Left/Right or Top/Bottom by comparing
        % the distance between x coords and y coords
        dx = abs(cx - ex);
        dy = abs(cy - ey);
        if(dy < dx)
            % grabbed a left/right edge point
            set(hp,'XData',cpt(1,1))
            [xarr,yarr] = get_current_edge_pts;
            minx = min(xarr);
            maxx = max(xarr);
            cx = 0.5*(minx + maxx);
            % update midpt location of edges top/bottom
            %    hp_edges are organized as L, R, B, T
            set(ud.hp_edges(3),'xdata',cx)
            set(ud.hp_edges(4),'xdata',cx)
        else
            % grabbed a top/bottom edge point
            set(hp,'YData',cpt(1,2))
            [xarr,yarr] = get_current_edge_pts;
            miny = min(yarr);
            maxy = max(yarr);
            cy = 0.5*(miny + maxy);
            % update midpt location of edges top/bottom
            set(ud.hp_edges(1),'ydata',cy)
            set(ud.hp_edges(2),'ydata',cy)
        end
        
        [xarr,yarr] = get_current_edge_pts;
        
        update_ctr_from_edges;
        update_corners_from_ctr(cx,cy);
        update_patch
        update_textbox_pos;
        
    end

    function move_ctr(varargin)
        oldcx = get(ud.hp_ctrpt,'xdata');
        oldcy = get(ud.hp_ctrpt,'ydata');
        cpt = get(ud.hax,'CurrentPoint');
        newcx = cpt(1,1);
        newcy = cpt(1,2);
        
        set(ud.hp_ctrpt,'xdata',newcx,'ydata',newcy)
        update_corners_from_ctr(newcx,newcy);
        update_edges_from_corners;
        update_patch;
        update_textbox_pos;
    end

    function move_corner(varargin)
        hp = gco;
        cpt = get(ud.hax,'CurrentPoint');
        % find which patch point is closest to this corner
        px = get(ud.hpatch,'xdata');
        py = get(ud.hpatch,'ydata');
        dist = sqrt((px - cpt(1,1)).^2 + (py - cpt(1,2)).^2);
        [junk,idx] = min(dist);
        % 
        newpx = px;
        newpy = py;
        % update only the point
        %                 BL BR UR UL
        %    px format:   x1 x2 x2 x1
        %    py format:   y1 y1 y2 y2
        %                 1  2  3  4
        if(idx == 1)
            newpx(1) = cpt(1,1);
            newpx(4) = cpt(1,1);
            newpy(1) = cpt(1,2);
            newpy(2) = cpt(1,2);
        elseif(idx == 2)
            newpx(2) = cpt(1,1);
            newpx(3) = cpt(1,1);
            newpy(1) = cpt(1,2);
            newpy(2) = cpt(1,2);
        elseif(idx == 3)
            newpx(3) = cpt(1,1);
            newpx(2) = cpt(1,1);
            newpy(3) = cpt(1,2);
            newpy(4) = cpt(1,2);
        elseif(idx == 4)
            newpx(4) = cpt(1,1);
            newpx(1) = cpt(1,1);
            newpy(4) = cpt(1,2);
            newpy(3) = cpt(1,2);
        end
        
        [newpx,newpw] = organize_corners(newpx,newpy);
        
        set(ud.hpatch,'xdata',newpx,'ydata',newpy)
        for k = 1:4
            set(ud.hp_corners(k),'xdata',newpx(k),'ydata',newpy(k))
        end
        
        % update all other objects
        update_edges_from_corners
        update_ctr_from_corners
        update_patch
        update_textbox_pos
    end

    function move_textbox(varargin)
        hax = get(ud.htext,'Parent');
        cpt = get(hax,'CurrentPoint');
        set(ud.htext,'Position',[cpt(1,1) cpt(1,2) 0])
        
        % update all other objects
        compute_textbox_offset_fig
        update_textbox_pos
    end

    function compute_textbox_offset_fig()
        [cx,cy] = get_current_ctr_pt;
        c_fig  = pointInAxesToFigurePosition([cx cy],ud.hax);
        cx_fig = c_fig(1);
        cy_fig = c_fig(2);
        tbpos = get(ud.htext,'Position');
        t_fig = pointInAxesToFigurePosition(tbpos(1:2), ud.hax);
        tx_fig = t_fig(1);
        ty_fig = t_fig(2);
        ud.textbox_offset_fig = [tx_fig-cx_fig  ty_fig-cy_fig];
    end

    function update_ctr_from_corners()
        [cx,cy] = compute_center_from_corners;
        set(ud.hp_ctrpt,'xdata',cx,'ydata',cy)
    end

    function update_ctr_from_edges()
        [xarr,yarr] = get_current_edge_pts;
        minx = min(xarr);
        maxx = max(xarr);
        miny = min(yarr);
        maxy = max(yarr);
        cx = 0.5*(minx + maxx);
        cy = 0.5*(miny + maxy);
        set(ud.hp_ctrpt,'xdata',cx,'ydata',cy);
    end

    function update_edges_from_corners()
        
        edges = compute_edges_from_corners;
        
        set(ud.hp_edges(1),'xdata',edges(1,1),'ydata',edges(1,2),'tag','left edge')
        set(ud.hp_edges(2),'xdata',edges(2,1),'ydata',edges(2,2),'tag','right edge')
        set(ud.hp_edges(3),'xdata',edges(3,1),'ydata',edges(3,2),'tag','bottom edge')
        set(ud.hp_edges(4),'xdata',edges(4,1),'ydata',edges(4,2),'tag','top edge')
    end

    function update_corners_from_ctr(ctrx,ctry)
        % called from move_ctr only!
        
        [xarr,yarr] = get_current_edge_pts;
        minx = min(xarr);
        maxx = max(xarr);
        miny = min(yarr);
        maxy = max(yarr);
        wx = abs(maxx - minx);
        wy = abs(maxy - miny);
        set(ud.hp_corners(1),'xdata',ctrx-wx/2,'ydata',ctry-wy/2)
        set(ud.hp_corners(2),'xdata',ctrx+wx/2,'ydata',ctry-wy/2)
        set(ud.hp_corners(3),'xdata',ctrx+wx/2,'ydata',ctry+wy/2)
        set(ud.hp_corners(4),'xdata',ctrx-wx/2,'ydata',ctry+wy/2)
    end

    function update_textbox_stats()
        set(ud.htext,'String',ud.stat_str)
    end

    function update_textbox_pos()
        [cx,cy]     = get_current_ctr_pt;
        c_fig = pointInAxesToFigurePosition([cx cy],gca);
        cx_fig = c_fig(1);
        cy_fig = c_fig(2);
        newc_fig = [ud.textbox_offset_fig(1)+cx_fig ud.textbox_offset_fig(2)+cy_fig];
        newc_axe = pointInFigureToAxesPosition(newc_fig,gca);
        
        set(ud.htext,'Position',[newc_axe(1) newc_axe(2)], ...
            'String',ud.stat_str);
        update_textbox_background
    end

    function update_textbox_background(varargin)
        tbpos = get(ud.htext,'Extent');
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
        
        set(ud.htextBG,'XData',tbBGposx,'YData',tbBGposy,'Visible','on')
    end

    function button_release(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',[])
        set(ud.hfig,'WindowButtonUpFcn',[])
        set(ud.hfig,'Interruptible','on')
    end

    function edges = compute_edges_from_corners()
        [xarr,yarr] = get_current_corner_pts;
        minx = min(xarr);
        maxx = max(xarr);
        miny = min(yarr);
        maxy = max(yarr);
        midx = 0.5*(maxx + minx);
        midy = 0.5*(maxy + miny);
        edges(1,1:2) = [minx midy]; % left edge
        edges(2,1:2) = [maxx midy]; % right edge
        edges(3,1:2) = [midx miny]; % bottom edge
        edges(4,1:2) = [midx maxy]; % top edge
    end

    function [cx,cy] = compute_center_from_corners()
        [xarr,yarr] = get_current_corner_pts;
        minx = min(xarr);
        maxx = max(xarr);
        miny = min(yarr);
        maxy = max(yarr);
        midx = 0.5*(maxx + minx);
        midy = 0.5*(maxy + miny);
        cx = midx;
        cy = midy;
    end
    
    function cleanup(varargin)
        set(ud.hpatch,'DeleteFcn',[])
        set(ud.htext,'DeleteFcn',[])
        set(ud.htextBG,'DeleteFcn',[])
        set(ud.hp_edges,'DeleteFcn',[])
        set(ud.hp_corners,'DeleteFcn',[])
        set(ud.hp_ctrpt,'DeleteFcn',[])
        delete(ud.hpatch)
        delete(ud.htext)
        delete(ud.htextBG)
        delete(ud.hp_edges)
        delete(ud.hp_corners)
        delete(ud.hp_ctrpt)
    end
    
    
end