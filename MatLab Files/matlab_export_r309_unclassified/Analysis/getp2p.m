function y = getp2p(varargin)

    if(nargin > 0)
        if(nargin == 1)
            xscale = varargin{1};
        elseif(nargin == 2)
            xscale = varargin{1};
            yscale = varargin{2};
        end
    else
        xscale = 1;
        yscale = 1;
    end
    u = ginput(2);
    dx = abs(diff(u(:,1)))*xscale;
    dy = diff(u(:,2))*yscale;
    m = dy/dx;
    msgbox(sprintf(' dx=%e\n 1/dx=%e\n dy=%e\n m = %e',dx,1/dx,dy,m),'Peak-to-Peak');
    
end