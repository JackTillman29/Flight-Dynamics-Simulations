function cb_drawslope_hz_midas(varargin)
    pick_pts = ginput(2);
    x1 = pick_pts(1,1);
    x2 = pick_pts(2,1);
    
    tMix = evalin('base','tMix');
    phaseRel = evalin('base','phaseRel');
    
    idx_data = find( (tMix >= x1) & (tMix <= x2) );
    y_data = phaseRel(idx_data);
    x_data = tMix(idx_data).';
    idx_n0 = find(~isnan(y_data));
    y_data = y_data(idx_n0);
    x_data = x_data(idx_n0);
    pres = polyfit(x_data,y_data,1);
    m = pres(1);
    b = pres(2);
    line_x = [x1 x2];
    line_y = [ ...
        m*x1 + b ...
        m*x2 + b];
    hLine = line(line_x,line_y);
    set(hLine,'Color','r','LineWidth',3);
    hText = text(x2,double(line_y(end)),sprintf(' %4.1f Hz',m/(2*pi)));
    set(hText,'Color','r','FontWeight','bold','BackgroundColor','w')
end

