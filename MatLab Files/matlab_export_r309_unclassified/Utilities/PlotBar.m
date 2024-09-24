function varargout = PlotBar(xcenter,xwidth,ymin,ymax,varargin)
    xw2 = xwidth / 2.0;
    h = patch([xcenter-xw2 xcenter-xw2 xcenter+xw2 xcenter+xw2],[ymin ymax ymax ymin],'k');
    %h = patch([xcenter-xw2 xcenter-xw2 xcenter+xw2 xcenter+xw2]+10,[ymin ymax ymax ymin]+.1,'k');
    for k = 1 : 2 : length(varargin)
        set(h,varargin{k},varargin{k+1});
    end
end