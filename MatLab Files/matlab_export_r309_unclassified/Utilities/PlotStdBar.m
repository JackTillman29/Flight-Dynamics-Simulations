function varargout = PlotStdBar(xmean,xstd,xthick,ypos,yheight,ythick,varargin)
    xw2 = xthick / 2.0;
    yw2 = ythick / 2.0;
    
    xp = [ ...
        xstd-xw2
        xstd-xw2
        xstd+xw2
        xstd+xw2
        xstd-xw2
        xstd-xw2
        0.5*xw2   % mean bar
        0.5*xw2
        -0.5*xw2
        -0.5*xw2   % end mean bar
        -xstd+xw2
        -xstd+xw2
        -xstd-xw2
        -xstd-xw2
        -xstd+xw2
        -xstd+xw2 
        -0.5*xw2   % mean bar
        -0.5*xw2
        0.5*xw2
        0.5*xw2   % end mean bar
        ];
    
    yp = [ ...
        yw2
        yheight
        yheight
        -yheight
        -yheight
        -yw2
        -yw2
        -0.5*yheight
        -0.5*yheight
        -yw2
        -yw2
        -yheight
        -yheight
        yheight
        yheight
        yw2 
        yw2
        0.5*yheight
        0.5*yheight
        yw2
        ];
    
    
    h = patch(xp+xmean,yp+ypos,'k');
    %h = patch([xcenter-xw2 xcenter-xw2 xcenter+xw2 xcenter+xw2]+10,[ymin ymax ymax ymin]+.1,'k');
    for k = 1 : 2 : length(varargin)
        set(h,varargin{k},varargin{k+1});
    end
    ht=text(min(xp+xmean),ypos-yheight,100,sprintf('%4.1f ', xmean-xstd));
    set(ht,'HorizontalAlignment','right');
    h2=text(xmean,ypos-yheight,100,sprintf('%4.1f', xmean));
    set(h2,'HorizontalAlignment','center');
    h3=text(xmean+xstd,ypos-yheight,100,sprintf('  %4.1f', xmean+xstd));
    set(h3,'HorizontalAlignment','left');
    
    set([ht h2 h3],'FontSize',7,'FontWeight','bold');
    
end