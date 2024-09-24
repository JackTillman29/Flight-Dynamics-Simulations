function add_secondary_units_to_tick_labels(haxes,select_axis,scale_function,labelstr)
% Syntax:
%   add_secondary_units_to_tick_labels(haxes,select_axis,scale_function,labelstr)
%       haxes:      axes handle
%       select_axis: 'x' or 'y' (single character)
%       scale_function: function handle of the math to do on the original
%                   axes data to arrive at the secondary unit string
%       labelstr:   character array of the units/label for the secondary
%                   tick label 
%
% author: Jeff Hole (Booz Allen Hamilton)
% 11-19-2019

    % if not already, set the TickLabelInterpreter to 'latex'
    set(haxes,'TickLabelInterpreter','latex');
    
    % for the callback update function (if using), save the inputs to the
    % UserData structure
    ud = get(haxes,'UserData');
    if(~isfield(ud,'ticklabel_secondary_units'))
        ud.ticklabel_secondary_units.select_axis = select_axis;
        ud.ticklabel_secondary_units.scale_function = scale_function;
        ud.ticklabel_secondary_units.labelstr = labelstr;
        set(haxes,'UserData',ud);
        
%         % now add an update callback context menu to the axes object
%         hcmenu = uicontextmenu();
%         uimenu(hcmenu, 'Label', 'Update secondary tick labels', ...
%             'Callback', @cb_add_secondary_units_to_tick_labels);
%         set(haxes,'UIContextMenu',hcmenu);
        
        % add a listener to the axes axis
        listenobj1 = addlistener(haxes,'XLim','PostSet',@cb_add_secondary_units_to_tick_labels);
        listenobj2 = addlistener(haxes,'YLim','PostSet',@cb_add_secondary_units_to_tick_labels);
    end
    
    drawnow;
    
    switch lower(select_axis)
        case 'x'
            ticks = get(haxes,'xtick');
        case 'y'
            ticks = get(haxes,'ytick');
    end
    
    for k = 1:length(xticks)
        ticklabs{k} = ['\begin{tabular}{c} ' num2str(ticks(k)) ' \\ (' num2str(scale_function(ticks(k))) ' ', labelstr ') \end{tabular}'];
    end
    
    switch lower(select_axis)
        case 'x'
            set(haxes,'xticklabels',ticklabs);
        case 'y'
            set(haxes,'yticklabels',ticklabs);
    end
    
    % this pause helps the callback feature use only the most up-to-date
    % information when computing tick label strings.
%     drawnow;
    pause(0.1)
    
    function varargout = cb_add_secondary_units_to_tick_labels(varargin)
    % haxes,select_axis,scale_factor,labelstr

    varargin{1}
    % varargin{2}

    % have to right-click on axes to get to this point, so we know gca is
    % correct axes object handle
    haxes = gca;

    ud = get(haxes,'UserData');
    if(isfield(ud,'ticklabel_secondary_units'))
        select_axis  = ud.ticklabel_secondary_units.select_axis;
        scale_function = ud.ticklabel_secondary_units.scale_function;
        labelstr     = ud.ticklabel_secondary_units.labelstr;
    end

    drawnow;
    add_secondary_units_to_tick_labels(haxes,select_axis,scale_function,labelstr)

    end
    
    
end
