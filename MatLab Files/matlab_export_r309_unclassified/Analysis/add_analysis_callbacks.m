function add_analysis_callbacks(varargin)
if(nargin > 0)
    hfig = varargin{1};
else
    hfig = gcf;
end
set(hfig,'Color','w');
% % Define a context menu; it is not attached to anything
% hcmenu = uicontextmenu;
% % Define callbacks for context menu items that change linestyle
% hcb1 = ['if(~exist(''pp'')) pp=PrepForPrint(); end; PrepForPrint(gcf,pp);'];
% % Define the context menu items and install their callbacks


%% define context menus
% if context menu already exists in figure, add to it
hcmenu = get(hfig,'uicontextmenu');
% figure doesn't contain a context menu object, create a new one.
if(isempty(hcmenu))
    hcmenu = uicontextmenu;
end
item1 = uimenu(hcmenu, 'Label', 'DrawSlope', 'Callback', 'cb_drawslope');
item2 = uimenu(hcmenu, 'Label', 'dx', 'Callback', 'cb_dx');
item3 = uimenu(hcmenu, 'Label', 'dy', 'Callback', 'cb_dy');
item4 = uimenu(hcmenu, 'Label', 'Save to workspace', 'Callback', 'cb_saveToWorkspace');
set(hfig,'UIContextMenu',hcmenu);



%% define tool bar
% find existing user-defined toolbar and add to it
htb = findobj('Parent',hfig,'Type','uitoolbar');
if(~isempty(htb))
    if(length(htb) > 1 & strcmp(get(hfig,'Toolbar'),'auto')) % create new toolbar (first one is the default toolbar)
        htb = uitoolbar(hfig);
        new_toolbar_created = 1;
    else % add to most recently created toolbar
        htb = htb(1);
        new_toolbar_created = 0;
    end
else
    new_toolbar_created = 1;
    htb = uitoolbar(hfig);
    % get directory where this function is stored (which should contain misc\)
    analysis_path = fileparts(mfilename('fullpath'));
    item1 = uipushtool(htb,'CData',imread([analysis_path,'/misc/dx_button.png']),...
        'TooltipString','Compute dx',...
        'ClickedCallback',@cb_dx);
    item1b = uipushtool(htb,'CData',imread([analysis_path,'/misc/frq_button.png']),...
        'TooltipString','Compute Frequency',...
        'ClickedCallback',@cb_frq);
    item2 = uipushtool(htb,'CData',imread([analysis_path,'/misc/dy_button.png']),...
        'TooltipString','Compute dy',...
        'ClickedCallback',@cb_dy);
    item3 = uipushtool(htb,'CData',imread([analysis_path,'/misc/m_button.png']),...
        'TooltipString','Compute slope',...
        'ClickedCallback',@cb_m);
    item4 = uipushtool(htb,'CData',imread([analysis_path,'/misc/move_data_2.png']),...
        'TooltipString','Move data to another figure',...
        'ClickedCallback',@cb_move);
    item5 = uipushtool(htb,'CData',imread([analysis_path,'/misc/move_data_waveObj.png']),...
        'TooltipString','Combine WaveObj data...',...
        'ClickedCallback',@cb_move_waveObj);
    item6 = uipushtool(htb,'CData',imread([analysis_path,'/misc/avg_window.png']),...
        'TooltipString','Compute Average Within Area (image)',...
        'ClickedCallback',@cb_two_point_average_box);
    item7 = uipushtool(htb,'CData',imread([analysis_path,'/misc/avg_region.png']),...
        'TooltipString','Compute Average Within Region (line)',...
        'ClickedCallback',@cb_two_point_average_line);
    item8 = uipushtool(htb,'CData',imread([analysis_path,'/misc/grab_slice_slant.png']),...
        'TooltipString','Grab Slice from Image Data',...
        'ClickedCallback',@cb_grab_slice_slant);
    item9 = uipushtool(htb,'CData',imread([analysis_path,'/misc/grab_slice_horizontal.png']),...
        'TooltipString','Grab Slice from Image Data',...
        'ClickedCallback',@cb_grab_slice_horizontal);
    item10 = uipushtool(htb,'CData',imread([analysis_path,'/misc/grab_slice_vertical.png']),...
        'TooltipString','Grab Slice from Image Data',...
        'ClickedCallback',@cb_grab_slice_vertical);
    item11 = uipushtool(htb,'CData',imread([analysis_path,'/misc/make_data_callout.png']),...
        'TooltipString','Make a data callout in figure',...
        'ClickedCallback',@cb_make_data_callout);
    % add fft tools to toolbar
    item12 = uipushtool(htb,'CData',imread([analysis_path,'/misc/fft_button.png']),...
        'TooltipString','Compute FFT of windowed data',...
        'ClickedCallback',@cb_fft);
end

%% define menu bar

% hmb = uimenu(hfig);



%% define keyboard callbacks
% THE FOLLOWING ONLY WORKS WITH LATER VERSIONS OF MATLAB (>2014)
set(hfig,'KeyPressFcn',@cb_modify_colorbar_caxis);

% set(hfig,'KeyPressFcn',@cb_keyPress);

%% finally, re-adjust figure position due to the offset created by adding toolbar
% (adding toolbar may push the figure up above the top of the screen, depending on figure defaults)
if(new_toolbar_created == 1)
    figpos = get(hfig,'Position');
    switch lower(get(hfig,'units'))
        case 'pixels'
            figpos(2) = figpos(2)-25;
        case 'normalized'
            figpos(2) = figpos(2)-0.02;
    end
    set(hfig,'Position',figpos)
end

end