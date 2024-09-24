function varargout = cb_plot_keypress(varargin)
    h = get(varargin{1},'UserData');
    switch varargin{2}.Key
        case 'uparrow'
            h.this.prevPage();
        case 'downarrow'
            h.this.nextPage();
        otherwise
    end
    plot(h.this,h.fcn,h.hi)
end