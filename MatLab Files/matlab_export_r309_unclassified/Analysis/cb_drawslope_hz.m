function cb_drawslope_hz(varargin)
    mode = 1;

    pick_pts = ginput(2);
    x1 = pick_pts(1,1);
    x2 = pick_pts(2,1);
    if(mode == 1)
        y1 = pick_pts(1,2);
        y2 = pick_pts(2,2);
    end
    
    
    if(mode == 1)
        x_data = [x1 x2];
        y_data = [y1 y2];
    else
        try
            t_on = evalin('base','t_on');
        catch
            t_on = evalin('caller','t_on')
        end
        refsig_on = evalin('base','refsig_on');

        idx_data = find( (1e6*t_on >= x1) & (1e6*t_on <= x2) );
        y_data = angle(refsig_on(idx_data));
        x_data = t_on(idx_data);
        idx_n0 = find(y_data ~= 0);
        
        y_data = y_data(idx_n0);
        x_data = x_data(idx_n0);
    end
    pres = polyfit(x_data,y_data,1);
    m = pres(1);
    b = pres(2);
    line_x = [x1 x2];
    if(mode == 1)
        line_y = [ ...
            m*x1 + b ...
            m*x2 + b];
    else
        line_y = [ ...
            1e-6*m*x1 + b ...
            1e-6*m*x2 + b];
    end
    hLine = line(line_x,line_y);
    set(hLine,'Color','b','LineWidth',3);
    if(mode == 1)
        hText = text(x2,double(line_y(end)),sprintf(' %4.1f Hz',1e6*m/(2*pi)));
    else
        hText = text(x2,double(line_y(end)),sprintf(' %4.1f Hz',m/(2*pi)));
    end
    set(hText,'Color','b','FontWeight','bold','BackgroundColor','w')
    colors = jet(30);
    for i = 1:3
        for j = 1:30
            set(hLine,'Color',colors(j,:))
            set(hText,'Color',colors(j,:))
            pause(0.01)
        end
    end
    set(hLine,'Color','b')
    set(hText,'Color','b')
    
end