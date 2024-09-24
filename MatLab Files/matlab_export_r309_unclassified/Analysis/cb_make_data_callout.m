function cb_make_data_callout(varargin)
% author: Jeff Hole (Booz Allen Hamilton)
% 12-1-2019

hfig = gcf;

h_pushtool = varargin{1};
h_actiondata = varargin{2};

% look for a default structure containing various properties used to build
% all data callout instances
does_exist = evalin('base','exist(''cp'')');
if(~does_exist)
    callout_props = PrepDataCallout;
else
    callout_props = evalin('base','cp');
end

% hax_main = gca

% get userdata from active figure (going to add a field tp this)
usrdat = get(hfig,'UserData');

% disp('click on the main axes you want to call out data')
% dontcare = ginput(1);
% 
% hax_main = gca;
% 

% return
% check if the following are on the MATLAB search path:
%   pointInFigureToAxesPosition
%   pointInAxesToFigurePosition
%   convertLBWHToPoints
%  ** (all should be in MATLAB\Utilities)

if(~exist('pointInFigureToAxesPosition','file'))
    res=which('cb_make_data_callout');
    try
        addpath(strrep(res,'\Analysis\cb_make_data_callout.m','\Utilities'));
    catch
        error('need to add MATLAB\Utilities to the path')
    end
end

%==========================================================================
% MAKE CALLOUT
%==========================================================================

% ENSURE THAT FIGURE UNITS ARE NORMALIZED (
original_figure_units_setting = get(hfig,'Units');
set(hfig,'units','normalized'); % make sure the figure units is set to 'normalized'

% SELECT ZOOM REGION - first use the zoom tool to select a zoom region
before_zoom_xlim = xlim;
before_zoom_ylim = ylim;
zoom(hfig)
disable_all_listeners(hfig)
waitfor(hfig,'CurrentPoint');
hax_main = gca; % get the current axes AFTER the click occurs
waitfor(hax_main,'xlim'); % wait for the xlim to change

% get handles in current axes that are to be copied to callout axes
hchildren_to_copy = get(hax_main,'children');

    % DISABLE ALL LISTENERS (IF THEY EXIST)
    disable_all_listeners(hfig)
after_zoom_xlim = xlim;
after_zoom_ylim = ylim;
xlim(hax_main,before_zoom_xlim);
ylim(hax_main,before_zoom_ylim);

% order: bottom-left, bottom-right, top-right, top-left
patx = [after_zoom_xlim(1) after_zoom_xlim(2) after_zoom_xlim(2) after_zoom_xlim(1)];
paty = [after_zoom_ylim(1) after_zoom_ylim(1) after_zoom_ylim(2) after_zoom_ylim(2)];

% FRAME (AXES -> AXES (NORMALIZED) COORDINATES)
% box size in format: [left bottom width height]
frame_pts_in_axe = [patx.' paty.'];
frame_pos_in_axe = [min(patx) min(paty) abs(diff(patx(1:2))) abs(diff(paty(2:3)))];
frame_pts_in_fig = pointInAxesToFigurePosition(frame_pts_in_axe, hax_main);
frame_pos_in_fig = [min(frame_pts_in_fig(:,1)) min(frame_pts_in_fig(:,2)) ...
                    abs(diff(frame_pts_in_fig(1:2,1))) abs(diff(frame_pts_in_fig(2:3,2)))];
    
drawnow;
% draw frame around data-of-interest
% hpatch_data_of_interest_frame = patch(patx,paty,'k','facecolor','none','edgecolor','k');
zoom(hfig,'off')

% RE-ENABLE ALL LISTENERS (IF THEY EXIST)
    drawnow;
    enable_all_listeners(hfig)

% NOW SELECT PLACE ON FIGURE TO DRAW NEXT AXES OBJ WITH ZOOMED IN PORTION
set(hfig,'Units','normalized')
k = waitforbuttonpress;
callout_pos_in_fig = rbbox;
callout_pts_in_fig = convertLBWHToPoints(callout_pos_in_fig);



%==========================================================================
%==========================================================================
%==========================================================================
% GET KEY VALUES OUT OF THE MAIN AXES OBJECT
%==========================================================================
%==========================================================================
%==========================================================================

% find if any image/surf plots exist in the figure
tmpchildlist = flipud(hchildren_to_copy); % flip to get oldest objects first

% look for the first image/surface object
for k = 1:length(hchildren_to_copy)
    hchild = hchildren_to_copy(k);
    switch get(hchild,'Type')
        case 'image'
            hdata = hchild;
            break
        case 'surface'
            hdata = hchild;
            break
    end
end
main_props_caxis = caxis(hax_main);
% main_props_alim = get(hax_main,'ALim');


% add a field to the UserData structure in hax_main
hax_main_ud = get(hax_main,'UserData');
hax_main_ud.data_callout = [];
hax_main_ud.data_callout.axes_object_label = 'main_axes';
set(hax_main,'UserData',hax_main_ud)



%==========================================================================
% DRAW SHADOW FIRST
%==========================================================================
% CALLOUT (FIGURE -> AXES COORDINATES)
callout_pts_axecoords = pointInFigureToAxesPosition(callout_pts_in_fig,hax_main);
frame_pts_axecoords   = convertLBWHToPoints(frame_pos_in_axe);
combined_pts_in_axe   = [callout_pts_axecoords; frame_pts_axecoords];
[K,V] = convhull( combined_pts_in_axe(:,1), combined_pts_in_axe(:,2) );

context_menu_obj_shadow = uicontextmenu;
m1 = uimenu(context_menu_obj_shadow,'Label','apply visual properties', ...
    'Callback',@cb_make_data_callout__update_properties);

hpatch_shadow = patch(hax_main,...
    combined_pts_in_axe(K,1),combined_pts_in_axe(K,2),callout_props.shadow_color, ...
    'Clipping','off', ...
    'FaceAlpha',callout_props.shadow_alpha,...
    'EdgeColor',callout_props.shadow_edge_color, ...
    'LineWidth',callout_props.shadow_edge_linewidth, ...
    'UIContextMenu',context_menu_obj_shadow);




%==========================================================================
% DRAW DATA FRAME
%==========================================================================
context_menu_obj_frame = uicontextmenu;
m1 = uimenu(context_menu_obj_frame,'Label','update',...
    'Callback',@cb_make_data_callout__update);
m2 = uimenu(context_menu_obj_frame,'Label','apply visual properties', ...
    'Callback',@cb_make_data_callout__update_properties);
m3 = uimenu(context_menu_obj_frame,'Label','make invisible', ...
    'Callback',@cb_make_data_callout__make_invisible);
m4 = uimenu(context_menu_obj_frame,'Label','make visible', ...
    'Callback',@cb_make_data_callout__make_visible);

hax_frame = axes('Position',frame_pos_in_fig, ...
    'Color',get(hax_main,'color'),...
    'xtick',[],'ytick',[],'box','on',...
    'xcolor',callout_props.axes_border_color, ...
    'ycolor',callout_props.axes_border_color, ...
    'linewidth',callout_props.axes_border_linewidth, ...
    'UIContextMenu',context_menu_obj_frame);

hax_frame_ud = get(hax_frame,'UserData');
hax_frame_ud.data_callout = [];
hax_frame_ud.data_callout.axes_object_label = 'frame_axes';
set(hax_frame,'UserData',hax_frame_ud)

copyobj(hchildren_to_copy,hax_frame);
xlim(hax_frame,frame_pos_in_axe(1) + [0 frame_pos_in_axe(3)])
ylim(hax_frame,frame_pos_in_axe(2) + [0 frame_pos_in_axe(4)])
caxis(hax_frame,main_props_caxis)





%==========================================================================
% FINALLY, DRAW CALLOUT
%==========================================================================
context_menu_obj_callout = uicontextmenu;
m1 = uimenu(context_menu_obj_callout,'Label','update', ...
    'Callback',@cb_make_data_callout__update);%_from_callout_limits);
m2 = uimenu(context_menu_obj_callout,'Label','apply visual properties', ...
    'Callback',@cb_make_data_callout__update_properties);

% giving this axes a string in "UserData" property to be used during
% interactive updates

hax_callout = axes('Position',callout_pos_in_fig, ...
    'Parent',hfig, ...
    'Color',get(hax_main,'color'), ...
    'xtick',[],'ytick',[],'box','on', ...
    'xcolor',callout_props.axes_border_color, ...
    'ycolor',callout_props.axes_border_color, ...
    'linewidth',callout_props.axes_border_linewidth, ...
    'UIContextMenu',context_menu_obj_callout);

hax_callout_ud = get(hax_callout,'UserData');
hax_callout_ud.data_callout = [];
hax_callout_ud.data_callout.axes_object_label = 'callout_axes';
set(hax_callout,'UserData',hax_callout_ud)

% copy original axes into new axes
copyobj(hchildren_to_copy,hax_callout)

% apply zoom limits to callout axes
xlim(hax_callout,after_zoom_xlim)
ylim(hax_callout,after_zoom_ylim)
caxis(hax_callout,main_props_caxis)

% PROBLEM WITH THIS IS THE INFINITE LOOP THESE LISTENERS WILL CAUSE
% addlistener(hax_main,'XLim','PostSet',@cb_make_data_callout__update)
% addlistener(hax_main,'YLim','PostSet',@cb_make_data_callout__update)



% save off handles
if(~isfield(usrdat,'data_callout_capability'))
    num_callouts = 1;
    icallout = 1;
else
    num_callouts = length(usrdat.data_callout_capability);
    icallout = num_callouts + 1;
end
usrdat.data_callout_capability(icallout).hax_main      = hax_main;
usrdat.data_callout_capability(icallout).hpatch_shadow = hpatch_shadow;
usrdat.data_callout_capability(icallout).hax_frame     = hax_frame;
usrdat.data_callout_capability(icallout).hax_callout   = hax_callout;

% add listeners (introduced in MATLAB 2008a)
usrdat.data_callout_capability(icallout).listen_hax_frame_pos    = addlistener(hax_frame,  'Position','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_callout_pos  = addlistener(hax_callout,'Position','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_frame_xlim   = addlistener(hax_frame,  'XLim','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_frame_ylim   = addlistener(hax_frame,  'YLim','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_callout_xlim = addlistener(hax_callout,'XLim','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_callout_ylim = addlistener(hax_callout,'YLim','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_main_xlim    = addlistener(hax_main,'XLim','PostSet',@cb_make_data_callout__update);
usrdat.data_callout_capability(icallout).listen_hax_main_ylim    = addlistener(hax_main,'YLim','PostSet',@cb_make_data_callout__update);

% usrdat = set_all_listeners_nonrecursive(usrdat,icallout);

% apply usrdat to "UserData"
set(hfig,'UserData',usrdat);

% SET FIGURE UNITS BACK TO WHAT IT WAS BEFORE CALLING THIS FUNCTION
set(hfig,'Units',original_figure_units_setting)





    function cb_make_data_callout__update(varargin)
        % updates the frame location

        % this callback is attached to the frame and callout axes
        if(nargin > 0)
            % not used, but making available.
            h_menu       = varargin{1};
            h_actiondata = varargin{2};
        end
        
%         % disable frame axes position listener while this callback
%         % operates on it.
%         disable_all_listeners(hfig)
        
%         disp('update called')

        % we know this because it was right-clicked for the context menu...
        hax = gca; % temporary axes handle, using "UserData" contained in figure for specific handles

        hax_ud = get(hax,'UserData');
        axes_object_label = hax_ud.data_callout.axes_object_label;
        if(strcmp(axes_object_label,'frame_axes'))
            % if clicked on frame axes, then update the frame
            update_mode = 'update self';
        elseif(strcmp(axes_object_label,'callout_axes'))
            % if clicked on callout axes, then update the frame with callout axes
            % limits
            update_mode = 'update frame using callout view';
        else
            update_mode = 'update frame and callout views';
        end

        hfig = get(hax,'Parent');

        usrdat = get(hfig,'UserData');

        % find which callout set is being modified
        icallout = make_data_callout__find_callout_set(hax,usrdat.data_callout_capability);
        if(icallout == 0)
            icallout = 1;
        end

        hpatch_shadow = usrdat.data_callout_capability(icallout).hpatch_shadow;
        hax_main      = usrdat.data_callout_capability(icallout).hax_main;
        hax_frame     = usrdat.data_callout_capability(icallout).hax_frame;
        hax_callout   = usrdat.data_callout_capability(icallout).hax_callout;

        callout_pos_in_fig = get(hax_callout,'Position'); % LBWH
        callout_pts_in_fig = convertLBWHToPoints(callout_pos_in_fig);
        callout_pts_in_axe = pointInFigureToAxesPosition(callout_pts_in_fig,hax_main);

        frame_pos_in_fig = get(hax_frame,'Position'); % LBWH
        frame_pts_in_fig = convertLBWHToPoints(frame_pos_in_fig);
        frame_pts_in_axe = pointInFigureToAxesPosition(frame_pts_in_fig,hax_main);
        frame_pos_in_axe = convertPointsToLBWH(frame_pts_in_axe); % LBWH

        % compute frame xy limits from frame axes position
        frame_xlims_from_axe_pos = frame_pos_in_axe(1) + [0 frame_pos_in_axe(3)];
        frame_ylims_from_axe_pos = frame_pos_in_axe(2) + [0 frame_pos_in_axe(4)];

        frame_xlims   = get(hax_frame,'xlim');
        frame_ylims   = get(hax_frame,'ylim');
        callout_xlims = get(hax_callout,'xlim');
        callout_ylims = get(hax_callout,'ylim');
        main_xlims    = get(hax_main,'xlim');
        main_ylims    = get(hax_main,'ylim');

        % if the axes position has changed, but the axes xy limits have
        % not, then the frame axes xy limits and the callout axes xy limits
        % should still be the same
        same_xlims = any(frame_xlims == callout_xlims);
        same_ylims = any(frame_ylims == callout_ylims);

        % disable frame axes position listener while this callback
        % operates on it.
        disable_all_listeners(hfig)

        switch update_mode
            case 'update self'
                if(same_xlims && same_ylims)
                    % IF ONLY THE FRAME AXES POSITION HAS CHANGED, THEN CALLOUT
                    % AXES AND FRAME AXES XY LIMITS SHOULD STILL BE THE SAME, SO
                    % UPDATE THE FRAME AXES XY LIMITS TO REFLECT THE NEW VIEW

                    % update frame limits
                    xlim(hax_frame, frame_xlims_from_axe_pos)
                    ylim(hax_frame, frame_ylims_from_axe_pos)

                    % update callout xy limits
                    xlim(hax_callout, frame_xlims_from_axe_pos)
                    ylim(hax_callout, frame_ylims_from_axe_pos)
                else
                    % FRAME AXES POSITION HAS NOT CHANGED, BUT THE XY LIMITS HAVE,
                    % SO UPDATE THE POSITION OF THE FRAME AXES SO THAT THE FRAME
                    % AND MAIN AXES DATA ALIGNS

                    % convert frame xy limits to a series of 4 points that draw a
                    % box
                    frame_limits_pts_in_axe(:,1) = [frame_xlims(1) frame_xlims(2) frame_xlims(2) frame_xlims(1)].';
                    frame_limits_pts_in_axe(:,2) = [frame_ylims(1) frame_ylims(1) frame_ylims(2) frame_ylims(2)].';
                    frame_limits_pts_in_fig = pointInAxesToFigurePosition(frame_limits_pts_in_axe,hax_main);
                    frame_limits_pos_in_fig = convertPointsToLBWH(frame_limits_pts_in_fig);

                    set(hax_frame,'Position',frame_limits_pos_in_fig)
                    % update callout xy limits
                    xlim(hax_callout, frame_xlims)
                    ylim(hax_callout, frame_ylims)
                end
                
                % TODO: check if callout is also a main axes and update all
                % children "data_callout"s
                
            case 'update frame using callout view'
                % need to update the frame axes location and its xy limits
                callout_limits_pts_in_axe(:,1) = [callout_xlims(1) callout_xlims(2) callout_xlims(2) callout_xlims(1)].';
                callout_limits_pts_in_axe(:,2) = [callout_ylims(1) callout_ylims(1) callout_ylims(2) callout_ylims(2)].';
                callout_limits_pts_in_fig = pointInAxesToFigurePosition(callout_limits_pts_in_axe,hax_main);
                callout_limits_pos_in_fig = convertPointsToLBWH(callout_limits_pts_in_fig);

                set(hax_frame,'Position',callout_limits_pos_in_fig)
                xlim(hax_frame, callout_xlims)
                ylim(hax_frame, callout_ylims)

        %         % update callout xy limits
        %         xlim(hax_callout, frame_xlims)
        %         ylim(hax_callout, frame_ylims)
            case 'update frame and callout views'
                % TODO: if the main axes is modified, update the callout
                % UPDATE FRAME FIRST THEN CALLOUT FRAME
                if(same_xlims && same_ylims)
                    % IF ONLY THE FRAME AXES POSITION HAS CHANGED, THEN CALLOUT
                    % AXES AND FRAME AXES XY LIMITS SHOULD STILL BE THE SAME, SO
                    % UPDATE THE FRAME AXES XY LIMITS TO REFLECT THE NEW VIEW

%                     % update frame limits
                    xlim(hax_frame, frame_xlims_from_axe_pos)
                    ylim(hax_frame, frame_ylims_from_axe_pos)
% 
%                     % update callout xy limits
                    xlim(hax_callout, frame_xlims_from_axe_pos)
                    ylim(hax_callout, frame_ylims_from_axe_pos)
                else
                    % FRAME AXES POSITION HAS NOT CHANGED, BUT THE XY LIMITS HAVE,
                    % SO UPDATE THE POSITION OF THE FRAME AXES SO THAT THE FRAME
                    % AND MAIN AXES DATA ALIGNS

                    % convert frame xy limits to a series of 4 points that draw a
                    % box
                    frame_limits_pts_in_axe(:,1) = [frame_xlims(1) frame_xlims(2) frame_xlims(2) frame_xlims(1)].';
                    frame_limits_pts_in_axe(:,2) = [frame_ylims(1) frame_ylims(1) frame_ylims(2) frame_ylims(2)].';
                    frame_limits_pts_in_fig = pointInAxesToFigurePosition(frame_limits_pts_in_axe,hax_main);
                    frame_limits_pos_in_fig = convertPointsToLBWH(frame_limits_pts_in_fig);

                    set(hax_frame,'Position',frame_limits_pos_in_fig)
                    % update callout xy limits
                    xlim(hax_callout, frame_xlims)
                    ylim(hax_callout, frame_ylims)

                end
%                 disp(['update 3-3-1'])
                % need to update the frame axes location and its xy limits
                callout_limits_pts_in_axe(:,1) = [callout_xlims(1) callout_xlims(2) callout_xlims(2) callout_xlims(1)].';
                callout_limits_pts_in_axe(:,2) = [callout_ylims(1) callout_ylims(1) callout_ylims(2) callout_ylims(2)].';
                callout_limits_pts_in_fig = pointInAxesToFigurePosition(callout_limits_pts_in_axe,hax_main);
                callout_limits_pos_in_fig = convertPointsToLBWH(callout_limits_pts_in_fig);

                set(hax_frame,'Position',callout_limits_pos_in_fig)
                xlim(hax_frame, callout_xlims)
                ylim(hax_frame, callout_ylims)
        end

        % re-enable frame axes position listener while this
        % callback operates on it.
        enable_all_listeners(hfig)
        
        % update shadow
        make_data_callout__update_shadow;
        
    end




    function cb_make_data_callout__update_properties(varargin)
        % this callback is attached to the callout axes.
        % this callback can also be called like a function (used this way in
        % PrepDataCallout.m)
        %     this function differentiates between the two by inspecting the number
        %     of arguments

        % 3 ways to use:
        %    1. as a callback (main usage)
        %    2. as a function passing in the figure handle
        %    3. as a function passing in the figure handle and the callout_props
        %           structure
        if(nargin > 0)
            if(nargin == 1)
                % 1 input which is a figure handle
                hfig = varargin{1};
                callout_props = evalin('base','cp');
            elseif(nargin == 2 && isstruct(varargin{2}))
                % 2 inputs, one a figure handle, the other a structure (not a
                % handle)
                hfig = varargin{1};
                callout_props = varargin{2};
            elseif(nargin == 2 && ~isstruct(varargin{2}))
                % 2 arguments, one a handle, the other NOT a structure (is a
                % matlab.ui.eventdata.ActionData), indicates this is a CALLBACK
                % call
                h_menu       = varargin{1};
                h_actiondata = varargin{2};

                hax = gca;
                % doesn't matter which axes, we just know if this callback is called, the
                % axes is in the figure-of-interest.
                hfig = get(hax,'Parent');

                % look for a default structure containing various properties used to build
                % all data callout instances
                does_exist = evalin('base','exist(''cp'')');
                if(~does_exist)
                    callout_props = PrepDataCallout;
                else
                    callout_props = evalin('base','cp');
                end
            end
        end

        usrdat = get(hfig,'UserData');

        % % find which callout set is being modified
        % icallout = make_data_callout__find_callout_set(hax,usrdat.data_callout_capability);

        for icallout = 1:length(usrdat.data_callout_capability)
        % UserData.data_callout_capability is created in "cb_make_data_callout"
        set(usrdat.data_callout_capability(icallout).hpatch_shadow, ...
            'FaceColor',callout_props.shadow_color, ...
            'FaceAlpha',callout_props.shadow_alpha, ...
            'EdgeColor',callout_props.shadow_edge_color, ...
            'LineWidth',callout_props.shadow_edge_linewidth)

        set(usrdat.data_callout_capability(icallout).hax_frame, ...
            'xcolor', callout_props.axes_border_color, ...
            'ycolor', callout_props.axes_border_color, ...
            'LineWidth', callout_props.axes_border_linewidth)

        set(usrdat.data_callout_capability(icallout).hax_callout, ...
            'xcolor', callout_props.axes_border_color, ...
            'ycolor', callout_props.axes_border_color, ...
            'LineWidth', callout_props.axes_border_linewidth)
        end


    end

    function make_data_callout__update_shadow(varargin)
        % this callback is attached to the callout axes
        if(nargin > 0)
            h_menu       = varargin{1};
            h_actiondata = varargin{2};
        end

        % we know this because it was right-clicked for the context menu...
        hax = gca;

        hfig = get(hax,'Parent');

        usrdat = get(hfig,'UserData');

        % find which callout set is being modified
        icallout = make_data_callout__find_callout_set(hax,usrdat.data_callout_capability);

        hpatch_shadow = usrdat.data_callout_capability(icallout).hpatch_shadow;
        hax_main      = usrdat.data_callout_capability(icallout).hax_main;
        hax_callout   = usrdat.data_callout_capability(icallout).hax_callout;
        hax_frame     = usrdat.data_callout_capability(icallout).hax_frame;

        callout_pos_in_fig = get(hax_callout,'Position'); % LBWH
        callout_pts_in_fig = convertLBWHToPoints(callout_pos_in_fig);
        callout_pts_in_axe = pointInFigureToAxesPosition(callout_pts_in_fig,hax_main);

        frame_pos_in_fig = get(hax_frame,'Position'); % LBWH
        frame_pts_in_fig = convertLBWHToPoints(frame_pos_in_fig);
        frame_pts_in_axe = pointInFigureToAxesPosition(frame_pts_in_fig,hax_main);

        % finally, update shadow
        combined_pts_in_axe = [callout_pts_in_axe; frame_pts_in_axe];
        [K,V] = convhull( combined_pts_in_axe(:,1), combined_pts_in_axe(:,2) );
        set(hpatch_shadow,'XData',combined_pts_in_axe(K,1),'YData',combined_pts_in_axe(K,2));

    end



    function iCalloutSet = make_data_callout__find_callout_set(hobj,data_callout_capability,is_main_axes)
        % one thing to keep in mind is that a callout can become the main axes for
        % another callout (callouts within callouts) so that duplicate handles will
        % be represented in the "data_callout_capability" structure array.
        %

        % cell array of strings, use this to search
        tmp = data_callout_capability(1);
        if(~exist('is_main_axes'))
            tmp = rmfield(data_callout_capability(1),'hax_main');
        end
        fielddat = fieldnames( tmp(1) );

        try
            % if MATLAB version does not support the syntax, then instantiate each
            % expected field.
            %     a.x    % ("x" is a field of structure "a")
            %     string_variable = 'x'
            %     a.(string_variable) ==evaluates to== a.x

            for k = 1:length(data_callout_capability)
                for j = 1:length(fielddat)
                    tmp = data_callout_capability(k).(fielddat{j});
                    if(tmp == hobj)
                        iCalloutSet = k;
                        return
                    end
                end
            end
        catch
            for k = 1:length(data_callout_capability)
                for j = 1:4
                    switch j
                        case 1
                            tmp = data_callout_capability(k).hax_frame;
                        case 2
                            tmp = data_callout_capability(k).hax_callout;
                        case 3
                            tmp = data_callout_capability(k).hpatch_shadow;
                    end

                    if(tmp == hobj)
                        iCalloutSet = k;
                        return
                    end
                end
            end

            % if get here, then didn't find the object "hobj" in
            % "data_callout_capability" structure
            error('did not find object in data structure');


        end
        
        % if get here, pass in a third arg (can be anything, only checks if it
        % exists) to compare input axes to all axes including main axes handle
        iCalloutSet = make_data_callout__find_callout_set(hobj,data_callout_capability,1);
        
    end

    function disable_all_listeners(hfig)
        ud = get(hfig,'UserData');
        if(isfield(ud,'data_callout_capability'))
            for icallout = 1:length(ud.data_callout_capability)
                % find all listener objects in each "data_callout_capability"
                % structure
                tmpstruct = ud.data_callout_capability(icallout);
                tmpfields = fieldnames(tmpstruct);
                for ifield = 1:length(tmpfields)
                    if(~isempty(regexpi(tmpfields{ifield},'^listen')))
                        usrdat.data_callout_capability(icallout).(tmpfields{ifield}).Enabled = 0;
                    end
                end
            end
        end
        set(hfig,'UserData',ud)
    end

    function enable_all_listeners(hfig)
        ud = get(hfig,'UserData');
        if(isfield(ud,'data_callout_capability'))
            for icallout = 1:length(ud.data_callout_capability)
                % find all listener objects in each "data_callout_capability"
                % structure
                tmpstruct = ud.data_callout_capability(icallout);
                tmpfields = fieldnames(tmpstruct);
                for ifield = 1:length(tmpfields)
                    if(~isempty(regexpi(tmpfields{ifield},'^listen')))
                        usrdat.data_callout_capability(icallout).(tmpfields{ifield}).Enabled = 1;
                    end
                end
            end
        end
        set(hfig,'UserData',ud)
    end

    function ud = set_all_listeners_nonrecursive(ud,icallout)
        if(isfield(ud,'data_callout_capability'))
            tmpstruct = ud.data_callout_capability(icallout);
            tmpfields = fieldnames(tmpstruct);
            for ifield = 1:length(tmpfields)
                if(~isempty(regexpi(tmpfields{ifield},'^listen')))
                    usrdat.data_callout_capability(icallout).(tmpfields{ifield}).Recursive = 0;
                end
            end
        end
    end

end
