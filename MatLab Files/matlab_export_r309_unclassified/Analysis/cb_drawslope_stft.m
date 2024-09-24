function cb_drawslope_stft(varargin)
    pick_pts = ginput(2);
    x1 = pick_pts(1,1);
    x2 = pick_pts(2,1);
    y1 = pick_pts(1,2);
    y2 = pick_pts(2,2);
    
    % get data in the axes
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');
    
    % compute slope
    dt = (y2 - y1);
    df = (x2 - x1);
    m = dt / df;
    
    % get time units from xlabel
    xlab = get(gca,'XLabel');
    ylab = get(gca,'YLabel');
    xlab = xlab.String;
    ylab = ylab.String;
    xidx1 = strfind(xlab,'(');
    yidx1 = strfind(ylab,'(');
    xidx2 = strfind(xlab,')');
    yidx2 = strfind(ylab,')');
    xlab = xlab(xidx1:xidx2);
    ylab = ylab(yidx1:yidx2);
    
    if(strcmpi(ylab,'s'))
        timeunits = 1;
    elseif(strcmpi(ylab,'ms'))
        timeunits = 1e-3;
    elseif(strcmpi(ylab,'\mus'))
        timeunits = 1e-6;
    elseif(strcmpi(ylab,'ns'))
        timeunits = 1e-9;
    end
    timestr = ylab;
    
    if(strcmpi(xlab,'Hz'))
        frequnits = 1;
    elseif(strcmpi(xlab,'kHz'))
        frequnits = 1e3;
    elseif(strcmpi(xlab,'MHz'))
        frequnits = 1e6;
    elseif(strcmpi(xlab,'GHz'))
        frequnits = 1e9;
    end
    freqstr = xlab;
    
    xline = [x1 x2];
    yline = [y1 y2];
    hLine = line(xline,yline);
    set(hLine,'Color','r','LineWidth',3)
    hText = text(x2,y2,sprintf([' %4.1f ',freqstr],df));
    set(hText,'Color','b','FontWeight','bold','BackgroundColor','w')
    
    
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
end