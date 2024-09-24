function cb_two_point_average_line(varargin)
    x=ginput(2);
    hdat=findobj('Type','line','Parent',gca);
    if(isempty(hdat))
        % if not a line, check for scatter plot
        hdat=findobj('Type','scatter','Parent',gca);
    end
    ud.hax = gca;
    ud.hfig = get(ud.hax,'Parent');
    
% %     xa=find(get(hdat,'XData') > x(1,1),1,'first');
% %     xb=find(get(hdat,'XData') < x(2,1),1,'last');
% %     
% %     xstep = sign(xb - xa);
% % 
% %     % ensure data are in row vectors
% %     xdat = reshape(get(hdat,'XData'),1,[]);
% %     c = reshape(get(hdat,'YData'),1,[]);
% %     ud.c_window = c(xa:xstep:xb);
% %     xdat_window = xdat(xa:xstep:xb);

    xstart = min(x(:,1));
    xstop  = max(x(:,1));
    fetch_data(xstart,xstop);



    
    ylims = ylim;

    % if linear
    %avg1 = mean(c_window(:));

    % if log
    %avg2 = 10*log10(mean(10.^(0.1*c_window(:))));
    %disp('lin,log');
    %disp([avg1,avg2])
    %msgbox(sprintf('Linear Average: %18.6e\nLog Average: %8.3f',avg1,avg2),'Average');
%     hp = patch( ...
%         [x(1,1) x(1,1) x(2,1) x(2,1)], ...
%         [x(1,2) x(2,2) x(2,2) x(1,2)], ...
%         'r','FaceAlpha',0.2);
    hpatch = patch( ...
        [x(1,1) x(1,1) x(2,1) x(2,1)], ...
        [ylims(1) ylims(2) ylims(2) ylims(1)], ...
        'r','FaceAlpha',0.2);
    
    set(hpatch,'EdgeColor','none');
    a=questdlg('Linear or Log (dB) Data?','Data Question','Linear','Log','Log');
    ud.data_kind = a;
    compute_stats;
    
    midx = mean(x(:,1));
    midy = mean(x(:,2));
    
    htext = text(midx,midy,ud.stat_str);
    set(htext,'Color','k','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
    
    % save off items for interactive callbacks
    ud.hdat = hdat;
    ud.htext = htext;
    ud.hpatch = hpatch;
    ud.selected = 0;
    
    hold on;
    ud.hp_edges(1) = plot(x(1,1),0.5*sum(ylims),'ksq','MarkerFaceColor','k');
    ud.hp_edges(2) = plot(x(2,1),0.5*sum(ylims),'ksq','MarkerFaceColor','k');
    ud.hp_ctrpt    = plot(0.5*(x(1,1)+x(2,1)), 0.5*sum(ylims),'ksq','MarkerFaceColor','k');
    hold off;
    
    % now that edges and ctrpt exist, call 
    compute_textbox_offset_fig
    
    % ADD INTERACTIVE ELEMENT (if MATLAB version supports necessary
    % functionality)
    % add buttondown callbacks and listener
    set(ud.hpatch,'DeleteFcn',@cleanup)
    set(ud.htext,'DeleteFcn',@cleanup)
    if exist('addlistener') && isprop(ud.hfig,'WindowButtonMotionFcn')
        set(ud.hpatch,'ButtonDownFcn',@click_patch)
        set(ud.hp_edges,'ButtonDownFcn',@click_edge)
        set(ud.hp_ctrpt,'ButtonDownFcn',@click_ctr)
        set(ud.htext,'ButtonDownFcn',@click_textbox)

        % add a listener to "ylim" to update the patch height
        ud.listen_ylim = addlistener(gca,'YLim','PostSet',@update_ylim);
        ud.listen_ylim.Recursive = 0;
    end
    
    
    function fetch_data(xstart,xstop)
        
        if(nargin == 0)
            [xarr,yarr] = get_current_edge_pts;
            xstart = xarr(1);
            xstop  = xarr(2);
        end
        
        xa=find(get(hdat,'XData') >= xstart,1,'first');
        xb=find(get(hdat,'XData') <= xstop ,1,'last');

        xstep = sign(xb - xa);

        % ensure data are in row vectors
        xdat = reshape(get(hdat,'XData'),1,[]);
        c = reshape(get(hdat,'YData'),1,[]);
        ud.c_window = c(xa:xstep:xb);
        ud.xdat_window = xdat(xa:xstep:xb);
    end
    
    function compute_stats()
        if(strcmpi(ud.data_kind,'Linear'))
            avgv = mean(ud.c_window(:));
            minv = min(ud.c_window(:));
            maxv = max(ud.c_window(:));
            dx = diff(ud.xdat_window);
            int_c = cumsum(ud.c_window(2:end)).*dx;
            min_int_c = min(int_c);
            max_int_c = max(int_c);
            fmt = 'min: %8.3e\nmax: %8.3e\navg: %8.3e';
            fmt = 'min: %8.3f\nmax: %8.3f\navg: %8.3f\nmin(integral): %8.3f\nmax(integral): %8.3f';
        else
            ud.c_window = 10.^(0.1*ud.c_window);
            avgv = 10*log10(mean(ud.c_window(:)));
            minv = 10*log10(min(ud.c_window(:)));
            maxv = 10*log10(max(ud.c_window(:)));
            dx = diff(ud.xdat_window);
            int_c = cumsum(ud.c_window(2:end)).*dx;
            min_int_c = 10*log10(min(int_c));
            max_int_c = 10*log10(max(int_c));
            fmt = 'min: %8.3f\nmax: %8.3f\navg: %8.3f';
            fmt = 'min: %8.3f\nmax: %8.3f\navg: %8.3f\nmin(integral): %8.3f\nmax(integral): %8.3f';
        end
        ud.stat_str = sprintf(fmt,minv,maxv,avgv,min_int_c,max_int_c);
    end

    function update_ylim(varargin)
        ylims = get(ud.hax,'ylim');
        my = 0.5*sum(ylims);
        set(ud.hp_edges(1),'ydata',my);
        set(ud.hp_edges(2),'ydata',my);
        set(ud.hp_ctrpt   ,'ydata',my);
        update_patch;
        update_textbox_pos;
    end

    function update_patch(varargin)
        [xarr,yarr] = get_current_edge_pts;
        ylims = get(ud.hax,'ylim');
        set(ud.hpatch,'xdata',[xarr(1) xarr(2) xarr(2) xarr(1)], ...
            'ydata',[ylims(1) ylims(1) ylims(2) ylims(2)])
    end

    function [xarr,yarr] = get_current_edge_pts()
        x1 = get(ud.hp_edges(1),'xdata');
        y1 = get(ud.hp_edges(1),'ydata');
        x2 = get(ud.hp_edges(2),'xdata');
        y2 = get(ud.hp_edges(2),'ydata');
        xarr = sort([x1 x2],'ascend');
        yarr = sort([y1 y2],'ascend');
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
        else
            ud.selected = 0;
            set(ud.hp_edges,'Visible','off');
            set(ud.hp_ctrpt,'Visible','off');
        end
    end

    function click_edge(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_edge)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
    end

    function click_ctr(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_ctr)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
    end

    function move_edge(varargin)
        hp = gco;
        cpt = get(ud.hax,'CurrentPoint');
        set(hp,'xdata',cpt(1,1))
        update_ctr
        update_patch
        fetch_data
        compute_stats
        update_textbox_stats
    end

    function move_ctr(varargin)
        hp = gco;
        cpt = get(ud.hax,'CurrentPoint');
        set(ud.hp_ctrpt,'xdata',cpt(1,1))
        update_edges
        update_patch
        fetch_data
        compute_stats
        update_textbox_stats
        update_textbox_pos
    end

    function update_ctr()
        [xarr,yarr] = get_current_edge_pts;
        mx = 0.5*sum(xarr); % midpoint between edges
        my = 0.5*sum(yarr); 
        set(ud.hp_ctrpt,'xdata',mx,'ydata',my)
        update_textbox_pos
    end

    function update_edges()
        [xarr,yarr] = get_current_edge_pts;
        [cx,cy] = get_current_ctr_pt;
        wx = abs(diff(xarr));
        wy = abs(diff(yarr));
        set(ud.hp_edges(1),'xdata',cx-wx/2,'ydata',cy)
        set(ud.hp_edges(2),'xdata',cx+wx/2,'ydata',cy)
    end

    function update_textbox_stats()
        set(ud.htext,'String',ud.stat_str)
    end

    function update_textbox_pos(varargin)
        [cx,cy]     = get_current_ctr_pt;
        c_fig = pointInAxesToFigurePosition([cx cy],gca);
        cx_fig = c_fig(1);
        cy_fig = c_fig(2);
        newc_fig = [ud.textbox_offset_fig(1)+cx_fig ud.textbox_offset_fig(2)+cy_fig];
        newc_axe = pointInFigureToAxesPosition(newc_fig,gca);
        
        set(ud.htext,'Position',[newc_axe(1) newc_axe(2)], ...
            'String',ud.stat_str);
    end


    function click_textbox(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',@move_textbox)
        set(ud.hfig,'WindowButtonUpFcn',@button_release)
    end

    function move_textbox(varargin)
        hax = get(ud.htext,'Parent');
        cpt = get(hax,'CurrentPoint');
        set(ud.htext,'Position',[cpt(1,1) cpt(1,2) 0])
        compute_textbox_offset_fig
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

    function button_release(varargin)
        set(ud.hfig,'WindowButtonMotionFcn',[])
        set(ud.hfig,'WindowButtonUpFcn',[])
    end

    function cleanup(varargin)
        set(ud.htext,'DeleteFcn',[])
        set(ud.hpatch,'DeleteFcn',[])
        set(ud.hp_edges,'DeleteFcn',[])
        set(ud.hp_ctrpt,'DeleteFcn',[])
        delete(ud.htext)
        delete(ud.hpatch)
        delete(ud.hp_edges)
        delete(ud.hp_ctrpt)
    end



end
