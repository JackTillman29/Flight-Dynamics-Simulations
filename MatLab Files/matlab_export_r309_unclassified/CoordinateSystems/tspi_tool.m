function tspi_tool()

% TO USE THE TSPI TOOL, NEED TO HAVE MATLAB R2013b OR LATER*
%   *R2013b introduced the "table" data type that is output from
%       readtable(), writetable()

close all; clear all; clc

screenSize = get(0,'ScreenSize');
screenSize = screenSize(3:4);

guiSize = [650 500];

h_start.hfig = figure( ...
    'Position',[15 screenSize(2)-guiSize(2)-100 guiSize], ...
    'Name','TSPI Analysis Tool', ...
    'NumberTitle','off','ToolBar','none','MenuBar','none','Resize','off', ...
    'DeleteFcn',@ui_cb__close_main_gui, ...
    'Tag','tspi_tool');

fig_background_color = get(h_start.hfig,'Color');

topbottom = 15;
leftright = 15;

defaults.ui_edit.infile = 'type in or browse for file';
defaults.ui_popup = 'load a file...';
defaults.ui_popup_units = {'deg','rad'};
defaults.tspi_tool_config_file_ext = 'ttconfig';
defaults.start_directory = pwd;

field_ht = 20;
text_width = guiSize(1)/8;
edit_width = guiSize(1)/4;
popup_units_width = 80;
infile_fontsize = 8;
infile_edit_width = guiSize(1)/2;
write_button_width = guiSize(1) / 4;
write_button_ht = 2*field_ht;
button_fontsize = 10;
obj_spacing = 20;

fontsize = 12;
horzalign = 'right';

kk_start = 1;
kk_step = 1.5;

kk_after_infile = 0;
kk_after_inputs = 0;      % init here to 0, gets overwritten after edit boxes are init
kk_after_inputs_step = 3; % distance between edit boxes and pushbuttons
kk_after_analyze_tspi_button = 0;
kk_after_analyze_tspi_button_step = 2.0;
kk_write_button_step = 2;

build_gui__infile_block();
build_gui__textboxes();
build_gui__editboxes();
build_gui__popup_colheaders();
build_gui__popup_units();

build_gui__analyze_tspi_button();
build_gui__write_buttons();

ui_util__disable_edit_boxes();
ui_util__disable_pushbuttons();

build_gui__menubar();
build_gui__draw_logos();

%======================================================================
% UI CONSTRUCTORS
%======================================================================
    function build_gui__menubar()
        m = uimenu(h_start.hfig,'Text','File');
        uimenu(m,'Text','Load Configuration File','Callback',@ui_cb_menubar__load_config);
        uimenu(m,'Text','Save Configuration File','Callback',@ui_cb_menubar__save_config);
    end
    function build_gui__infile_block()
        left_pos = leftright;
        obj_width = text_width;
        kk = kk_start;
        h_start.ui_text.infile = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Infile:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        left_pos = left_pos+obj_width+obj_spacing;
        h_start.ui_edit.infile = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht infile_edit_width field_ht], ...
            'String',defaults.ui_edit.infile,'HorizontalAlignment',horzalign,...
            'FontSize',infile_fontsize);
        left_pos = left_pos+infile_edit_width;
        browse_button_width = obj_width-10;
        h_start.ui_button_browse = uicontrol('Style','pushbutton','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht browse_button_width field_ht], ...
            'String','Browse...','HorizontalAlignment',horzalign,...
            'FontSize',button_fontsize,'Callback',@ui_cb_button_browse_infile);
        
        left_pos = left_pos+browse_button_width;
        browse_button_width = obj_width-20;
        h_start.ui_button_load = uicontrol('Style','pushbutton','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht browse_button_width field_ht], ...
            'String','Load','HorizontalAlignment',horzalign,...
            'FontSize',button_fontsize,'Callback',@ui_cb_button_load_infile);
        kk_after_infile = kk + kk_step;
    end
    function build_gui__analyze_tspi_button()
        kk = kk_after_inputs + kk_after_inputs_step;
        left_pos = 2*leftright;
        left_pos = guiSize(1)/2 - 0.5*write_button_width;
        h_start.ui_dependent_buttons.analyze_tspi = ...
            uicontrol('Style','pushbutton','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht write_button_width write_button_ht], ...
            'String','Analyze TSPI File','HorizontalAlignment','center',...
            'FontSize',button_fontsize, ...
            'Callback',@ui_cb_button_analyze_tspi);
        kk_after_analyze_tspi_button = kk;
    end
    function build_gui__textboxes()
        % BUILD TEXT BOXES
        left_pos = leftright;
        obj_width = text_width;
        kk = kk_after_infile + 1;
        h_start.ui_text.lat = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Lat:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.lon = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Lon:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.alt = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Alt:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.yaw = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Yaw:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.pitch = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Pitch:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.roll = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Roll:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.time = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Time:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
        kk = kk + kk_step;
        h_start.ui_text.run = uicontrol('Style','text','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','Run:','HorizontalAlignment',horzalign,...
            'FontSize',fontsize);
    end
    function build_gui__editboxes()
        % BUILD EDIT BOXES
        left_pos = leftright+text_width+obj_spacing;
        obj_width = edit_width;
        kk = kk_after_infile + 1;
        h_start.ui_edit.lat = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.lon = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.alt = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.yaw = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.pitch = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.roll = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.time = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        kk = kk + kk_step;
        h_start.ui_edit.run = uicontrol('Style','edit','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width field_ht], ...
            'String','','HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',@ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons);
        
        kk_after_inputs = kk;
    end
    function build_gui__popup_colheaders()
        % BUILD EDIT BOXES
        left_pos = leftright+text_width+obj_spacing+edit_width+obj_spacing;
        obj_width = edit_width;
        popup_ht = 25;
        kk = kk_after_infile + 1;
        h_start.ui_popup.lat = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'lat'});
        kk = kk + kk_step;
        h_start.ui_popup.lon = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'lon'});
        kk = kk + kk_step;
        h_start.ui_popup.alt = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'alt'});
        kk = kk + kk_step;
        h_start.ui_popup.yaw = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'yaw'});
        kk = kk + kk_step;
        h_start.ui_popup.pitch = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'pitch'});
        kk = kk + kk_step;
        h_start.ui_popup.roll = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'roll'});
        kk = kk + kk_step;
        h_start.ui_popup.time = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'time'});
        kk = kk + kk_step;
        h_start.ui_popup.run = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize,'Callback',{@ui_cb_popup,'run'});
    end
    function build_gui__popup_units()
        % BUILD EDIT BOXES
        left_pos = leftright+text_width+obj_spacing+edit_width+obj_spacing+edit_width+obj_spacing;
        obj_width = popup_units_width;
        popup_ht = 25;
        
        kk = kk_after_infile + 1;
        %         h_start.ui_popup_units.lat = uicontrol('Style','popup','Parent',h_start.hfig, ...
        %             'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
        %             'String',defaults.ui_popup_units,'HorizontalAlignment',horzalign,...
        %             'FontSize',fontsize,'Callback',{@ui_cb_popup_units,'lat'});
        kk = kk + kk_step;
        %         h_start.ui_popup_units.lon = uicontrol('Style','popup','Parent',h_start.hfig, ...
        %             'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
        %             'String',defaults.ui_popup_units,'HorizontalAlignment',horzalign,...
        %             'FontSize',fontsize,'Callback',{@ui_cb_popup_units,'lon'});
        kk = kk + kk_step;
        %         h_start.ui_popup_units.alt = uicontrol('Style','popup','Parent',h_start.hfig, ...
        %             'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
        %             'String',defaults.ui_popup_units,'HorizontalAlignment',horzalign,...
        %             'FontSize',fontsize,'Callback',{@ui_cb_popup_units,'alt'});
        kk = kk + kk_step;
        h_start.ui_popup_units.yaw = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup_units,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize);%,'Callback',{@ui_cb_popup_units,'yaw'});
        kk = kk + kk_step;
        h_start.ui_popup_units.pitch = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup_units,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize);%,'Callback',{@ui_cb_popup_units,'pitch'});
        kk = kk + kk_step;
        h_start.ui_popup_units.roll = uicontrol('Style','popup','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width popup_ht], ...
            'String',defaults.ui_popup_units,'HorizontalAlignment',horzalign,...
            'FontSize',fontsize);%,'Callback',{@ui_cb_popup_units,'roll'});
        
    end
    function build_gui__write_buttons()
        %         left_pos = 2*leftright;
        left_pos = guiSize(1)/2 - 0.5*write_button_width;
        obj_width = write_button_width;
        kk = kk_after_analyze_tspi_button + kk_after_analyze_tspi_button_step;
        h_start.ui_dependent_buttons.analyze_write_tspi = uicontrol('Style','pushbutton','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width write_button_ht], ...
            'String','Write TSPI','HorizontalAlignment','center',...
            'FontSize',button_fontsize,'Callback',@ui_cb_button_write_tspi);
        kk = kk + kk_write_button_step;
        h_start.ui_dependent_buttons.write_esams = uicontrol('Style','pushbutton','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width write_button_ht], ...
            'String','Write ESAMS','HorizontalAlignment','center',...
            'FontSize',button_fontsize,'Callback',@ui_cb_button_write_esams);
        kk = kk + kk_write_button_step;
        h_start.ui_dependent_buttons.write_suppressor= uicontrol('Style','pushbutton','Parent',h_start.hfig, ...
            'Position',[left_pos  guiSize(2)-topbottom-kk*field_ht obj_width write_button_ht], ...
            'String','Write Suppressor','HorizontalAlignment','center',...
            'FontSize',button_fontsize,'Callback',@ui_cb_button_write_suppressor);
    end
    function build_gui__draw_logos()
        % bottom left, bottom right
        imsize = [write_button_width write_button_width];
        left_picture_pos  = [3*leftright topbottom imsize];
        right_picture_pos = [guiSize(1)-3*leftright-imsize(1) topbottom imsize];
        
        h_start.ui_picture.project_axes = axes( ...
            'Units','pixels', ...
            'Parent',h_start.hfig, ...
            'Position',left_picture_pos);
        [A,map,transparency] = imread('logo_gnome.png');
        h_start.ui_picture.project_image = image(A, ...
            'Parent',h_start.ui_picture.project_axes, ...
            'AlphaData',transparency, ...
            'ButtonDownFcn',@ui_cb_clicked_img);
        axis(h_start.ui_picture.project_axes,'equal')
        set(h_start.ui_picture.project_axes, ...
            'XTick',[],'YTick',[],'Box','off', ...
            'XColor','none','YColor','none','color','none')
        
        h_start.ui_picture.branch_axes = axes( ...
            'Units','pixels', ...
            'Parent',h_start.hfig, ...
            'Position',right_picture_pos);
        [A,map,transparency] = imread('logo_aflcmc_ezja.png');
        h_start.ui_picture.branch_image = image(A, ...
            'Parent',h_start.ui_picture.branch_axes, ...
            'AlphaData',transparency, ...
            'ButtonDownFcn',@ui_cb_clicked_img);
        axis(h_start.ui_picture.project_axes,'equal')
        set(h_start.ui_picture.branch_axes, ...
            'XTick',[],'YTick',[],'Box','off', ...
            'XColor','none','YColor','none','color','none')
    end

    function ui_cb_clicked_img(varargin)
        him = varargin{1};
        hit = varargin{2};
        cdat = get(him,'cdata');
        cdat = mod(uint8(255) - cdat,uint8(255));
        set(him,'cdata',cdat)
    end

%======================================================================
% UI CALLBACKS
%======================================================================
    function ui_cb_menubar__save_config(varargin)
        % get a default string ready for the gui
        try
            % if a TSPI file has already been loaded, use the path and
            % filename as the default
            default_filename = h_start.tspiObj.filename;
            [filepath,filename,dontcare] = fileparts(default_filename);
            default_filename = [filepath '\' filename '.' defaults.tspi_tool_config_file_ext];
        catch
            default_filename = ['Untitled.' defaults.tspi_tool_config_file_ext];
        end
        
        % user-interactive file selection
        [filename,filepath,filtind] = uiputfile( ...
            default_filename,'Save TSPI Tool Config File');
        % handle "cancel" buttonpress
        if ~ischar(filename)
            warndlg('did not save the config file!')
            return
        end

        % write a text file containing the following:
        % save LLA, YPR, Time and Run editbox fields, and YPR units fields
        fid = fopen([filepath filename],'w+');
        
        fprintf(fid,'%s, %s\n', 'wgs84_lat', get(h_start.ui_edit.lat,'String'));
        fprintf(fid,'%s, %s\n', 'wgs84_lon', get(h_start.ui_edit.lon,'String'));
        fprintf(fid,'%s, %s\n', 'wgs84_alt', get(h_start.ui_edit.alt,'String'));
        fprintf(fid,'%s, %s\n', 'yaw',       get(h_start.ui_edit.yaw,'String'));
        fprintf(fid,'%s, %s\n', 'pitch',     get(h_start.ui_edit.pitch,'String'));
        fprintf(fid,'%s, %s\n', 'roll',      get(h_start.ui_edit.roll,'String'));
        fprintf(fid,'%s, %s\n', 'time',      get(h_start.ui_edit.time,'String'));
        fprintf(fid,'%s, %s\n', 'run',       get(h_start.ui_edit.run,'String'));
        
        fprintf(fid,'%s, %s\n', 'units_yaw',   defaults.ui_popup_units{get(h_start.ui_popup_units.yaw,'Value')});
        fprintf(fid,'%s, %s\n', 'units_pitch', defaults.ui_popup_units{get(h_start.ui_popup_units.pitch,'Value')});
        fprintf(fid,'%s, %s\n', 'units_roll',  defaults.ui_popup_units{get(h_start.ui_popup_units.roll,'Value')});
        fclose(fid);
    end
    function ui_cb_menubar__load_config(varargin)
        % user-interactive file selection
        [filename,filepath,filtind] = uigetfile( ...
            [defaults.start_directory '\*.' defaults.tspi_tool_config_file_ext],'MultiSelect','off');
        % handle "cancel" buttonpress
        if isnumeric(filename) & isnumeric(filepath)
            return
        else
            selected_file = [filepath filename];
        end
        
        fid = fopen([filepath filename],'r');
        
        % given an input string, find the index into the cell array
        % "defaults.ui_popup_units" to set the "Value" property in the
        % popup menu
        popup_units_handler = @(thisstr) ...
            find(not(cellfun('isempty',strfind(defaults.ui_popup_units,thisstr))));
        
        % apply config file to appropriate fields
        while(~feof(fid))
            % get line
            thisline = fgetl(fid);
            
            params = strsplit(thisline,',');
            params{2} = strtrim(params{2}); % remove leading/trailing spaces
            % parse line
            switch lower(params{1})
                case 'wgs84_lat'
                    set(h_start.ui_edit.lat,'String',params{2});
                case 'wgs84_lon'
                    set(h_start.ui_edit.lon,'String',params{2});
                case 'wgs84_alt'
                    set(h_start.ui_edit.alt,'String',params{2});
                case 'yaw'
                    set(h_start.ui_edit.yaw,'String',params{2});
                case 'pitch'
                    set(h_start.ui_edit.pitch,'String',params{2});
                case 'roll'
                    set(h_start.ui_edit.roll,'String',params{2});
                case 'time'
                    set(h_start.ui_edit.time,'String',params{2});
                case 'run'
                    set(h_start.ui_edit.run,'String',params{2});
                case 'units_yaw'
                    index = popup_units_handler(params{2});
                    if(index > 0)
                        set(h_start.ui_popup_units.yaw,'Value',index);
                    else
                        error(['units_yaw: "', params{2},'" is not a valid input'])
                    end
                case 'units_pitch'
                    index = popup_units_handler(params{2});
                    if(index > 0)
                        set(h_start.ui_popup_units.pitch,'Value',index);
                    else
                        error(['units_pitch: "', params{2},'" is not a valid input'])
                    end
                case 'units_roll'
                    index = popup_units_handler(params{2});
                    if(index > 0)
                        set(h_start.ui_popup_units.roll,'Value',index);
                    else
                        error(['units_roll: "', params{2},'" is not a valid input'])
                    end
                otherwise
                    error(['problem with config file ' selected_file '. "' params{1} '" is an unrecognized field'])
            end
        end
        % if we get here, we should have been successful loading the file
        ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons;
        
    end
    function ui_cb_popup(varargin)
        huiobj = varargin{1};
        actionData = varargin{2};
        popup_id_str = varargin{3};
        % apply the
        colhdrs = get(huiobj,'String');
        % if only a single column header found (??)
        if(ischar(colhdrs))
            if strcmpi(colhdrs,defaults.ui_popup)
                % one value in the popup menu is actually the default
                % string
                return
            else
                % only one column header found...
                set(h_start.ui_edit.(popup_id_str),'String',colhdrs)
            end
        else
            if ~strcmpi(colhdrs{1},defaults.ui_popup)
                put_in_editbox = colhdrs{get(huiobj,'Value')};
                % if LLA or YPR, then strip off the last character if it is
                % a number
                if ~isempty(str2num(put_in_editbox(end)))
                    put_in_editbox = put_in_editbox(1:(end-1));
                end
                set(h_start.ui_edit.(popup_id_str),'String',put_in_editbox)
            end
        end
        % check if all edit boxes are filled, then enable pushbuttons if so
        ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons();
    end

    function ui_cb_button_browse_infile(varargin)
        % user-interactive file selection
        [filename,filepath,filtind] = uigetfile([defaults.start_directory '\*'],'MultiSelect','off');
        % handle "cancel" buttonpress
        if isnumeric(filename) & isnumeric(filepath)
            return
        else
            selected_file = [filepath filename];
        end
        
        % save off this new directory
        defaults.start_directory = filepath;
        
        % set infile edit box to this file
        set(h_start.ui_edit.infile,'String',selected_file);
        
        % load the data (call the load button callback to avoid duplication
        % of code)
        ui_cb_button_load_infile();
    end
    function ui_cb_button_load_infile(varargin)
        % may not need if automatically loading in after browsing
        % get string from edit_infile box
        selected_file = get(h_start.ui_edit.infile,'String');
        if ~strcmpi(defaults.ui_edit.infile,selected_file)
            colhead = ui_util__get_table_colheaders(selected_file);
            ui_util__set_popup_colheaders(colhead);
            % enable the edit buttons
            ui_util__enable_edit_boxes
            % check to see if edit boxes are full before disabling the
            % pushbuttons (found this bug upon loading new files in a
            % single session of the tool)
            ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons()
        end
        
        if ~isfield(h_start,'tspiObj')
            if ~strcmp(selected_file,defaults.ui_edit.infile)
                h_start.tspiObj = tspi_class(selected_file);
            end
        else
            % see if this is the same file
            if strcmp(h_start.tspiObj.filename, selected_file)
                % do nothing, continue on to re-open tool
            else
                % not the same file, load the new one (after being asked
                % if you are done with the previous file)
                answer = questdlg('Are you done with previous file?','Load new TSPI file...', ...
                    'Yes','No','No');
                switch answer
                    case 'Yes'
                        h_start.tspiObj = tspi_class(selected_file);
                    case 'No'
                        % do nothing, continue on to re-open tool
                        % ensure the filename remains what it was
                        set(h_start.ui_edit.infile,'String',h_start.tspiObj.filename)
                    otherwise
                        set(h_start.ui_edit.infile,'String',h_start.tspiObj.filename)
                end
            end
        end
    end

    function ui_cb_button_analyze_tspi(varargin)
        
        % apply search strings fetched from the gui
        h_start.tspiObj.searchString_wgs84_lat = get(h_start.ui_edit.lat,'String');
        h_start.tspiObj.searchString_wgs84_lon = get(h_start.ui_edit.lon,'String');
        h_start.tspiObj.searchString_wgs84_alt = get(h_start.ui_edit.alt,'String');
        h_start.tspiObj.searchString_yaw       = get(h_start.ui_edit.yaw,'String');
        h_start.tspiObj.searchString_pitch     = get(h_start.ui_edit.pitch,'String');
        h_start.tspiObj.searchString_roll      = get(h_start.ui_edit.roll,'String');
        h_start.tspiObj.searchString_time      = get(h_start.ui_edit.time,'String');
        h_start.tspiObj.searchString_run       = get(h_start.ui_edit.run,'String');
        
        % these are not a part of the tool right now (always in deg)
        %h_start.tspiObj.units_lat = get(h_start.ui_popup_units.lat,'String');
        %h_start.tspiObj.units_lon = get(h_start.ui_popup_units.lon,'String');
        %h_start.tspiObj.units_alt = get(h_start.ui_popup_units.alt,'String');
        h_start.tspiObj.units_yaw   = defaults.ui_popup_units{get(h_start.ui_popup_units.yaw,'Value')};
        h_start.tspiObj.units_pitch = defaults.ui_popup_units{get(h_start.ui_popup_units.pitch,'Value')};
        h_start.tspiObj.units_roll  = defaults.ui_popup_units{get(h_start.ui_popup_units.roll,'Value')};
        
        numProblemRegions = fix_TSPI(h_start.tspiObj,'preprocess');
        %         numProblemRegions = fix_TSPI__fast_count_problem_regions();
        fix_TSPI(h_start.tspiObj,'platform',1,'run',1);
        
        if(any(numProblemRegions > 0))
            ui_util__disable_pushbuttons__write;
        else
            ui_util__enable_pushbuttons__write;
        end
        
    end

    function [plat_array,run_array] = ui_util__prep_for_write()
        nPlats = h_start.tspiObj.nPlats;
        nRuns = h_start.tspiObj.nRuns;
        
        % function ui_util__dialog_platform_run_selection()
        table_size = [400 400];
        screen_size = get(0,'ScreenSize'); screen_size = screen_size(3:4);
        screen_ctr = screen_size/2;
        leftright = 10;
        top = 10;
        bottom = 40;
        fig_size = table_size + [2*leftright top+bottom];
        fig_pos = [screen_ctr-fig_size/2 fig_size];
        hdialog = dialog('Units','pixels','Position',fig_pos);
        % put a matrix of logicals into a ColumnEditable uitable creates an
        % editable checkbox array!
        initial_data = logical(zeros(nRuns,nPlats));
        htab = uitable('Parent',hdialog,'Units','pixels', ...
            'Position',[leftright bottom table_size], ...
            'Data',initial_data, ...
            'ColumnEditable',true);
        for krun = 1:nRuns
            table_rownames{krun} = ['Run ',num2str(krun)];
        end
        for kplat = 1:nPlats
            table_colnames{kplat} = ['Plat ',num2str(kplat)];
        end
        set(htab,'RowName',table_rownames);
        set(htab,'ColumnName',table_colnames);
        set(htab,'ColumnWidth',num2cell(50*ones(1,nPlats)));
        
        button_ht = 30;
        button_bottom = (bottom-button_ht)/2;
        button_width = table_size(2)/2;
        button_pos = [fig_size(2)/2-button_width/2 button_bottom button_width button_ht];
        hbutton = uicontrol('Parent',hdialog, 'Position', button_pos, ...
            'String','OK','Callback',@cb_button_ok);
        
        output = nan;
        
        uiwait(hdialog)
        
        if(any(~isnan(output)))
            [run_array,plat_array] = find(output);
        else
            run_array = [];
            plat_array = [];
        end
        function cb_button_ok(varargin)
            % set output to the table
            output = get(htab,'Data');
            delete(hdialog)
        end
    end

    function ui_cb_button_write_esams(varargin)
        [plat_array,run_array] = ui_util__prep_for_write();
        % put files in the same place as the input file
        [path,filename,ext] = fileparts(h_start.tspiObj.filename);
        % call tspi_class method to output to ECEF for ESAMS
        if(isempty(plat_array) || isempty(run_array))
            return;
        end
        
        files_written = h_start.tspiObj.dumpESAMS_ECEF('platform',plat_array,'run',run_array);
        % make message box strings
        msgbox_cells{1} = 'Files written:';
        for k = 1:length(files_written)
            msgbox_cells{k+1} = files_written{k};
        end
        msgbox(msgbox_cells.')
    end



    function ui_cb_button_write_tspi(varargin)
        numProblemRegions = fix_TSPI(h_start.tspiObj,'preprocess');
        
        if(any(numProblemRegions > 0))
            answer = questdlg( ...
                {'Problems detected in TSPI file.','Still want to write a new TSPI file?','yeah that works'},'', ...
                'Yes','No','No');
            if ~strcmp(answer,'Yes')
                return
            end
        end
        
        % get a default string ready for the gui
        try
            % if a TSPI file has already been loaded, use the path and
            % filename as the default
            default_filename = h_start.tspiObj.filename;
            [filepath,filename,fileext] = fileparts(default_filename);
            default_filename = [filepath '\' filename '_fixed' fileext];
        catch
            default_filename = ['Untitled.' defaults.tspi_tool_config_file_ext];
        end
        
        % user-interactive file selection
        [filename,filepath,filtind] = uiputfile( ...
            default_filename,'Save TSPI File As');
        % handle "cancel" buttonpress
        if ~ischar(filename)
            warndlg('did not save the TSPI file!')
            return
        end
        writetable(h_start.tspiObj.tableData,[filepath filename])
    end
    function ui_cb_button_write_suppressor(varargin)
        disp('ui_cb_button_write_suppressor')
    end

    function ui_cb_edits__check_all_edits_and_update_analyze_pushbuttons(varargin)
        fields = fieldnames(h_start.ui_edit);
        populated = zeros(1,length(fields)+1);
        for k = 1:length(fields)
            % check if something is in the edit box
            populated(k) = ~isempty(get(h_start.ui_edit.(fields{k}),'String'));
        end
        populated(k+1) = exist(get(h_start.ui_edit.infile,'String'),'file');
        % enable the analyze pushbutton
        if all(populated)
            ui_util__enable_pushbuttons__analyze;
        else
            ui_util__disable_pushbuttons__analyze;
        end
    end

%======================================================================
% UI UTILITIES
%======================================================================
    function column_headers = ui_util__get_table_colheaders(file)
        tmptbl = readtable(file);
        column_headers = tmptbl.Properties.VariableNames;
    end
    function ui_util__set_popup_colheaders(column_headers)
        % column_headers needs to be a cell array of strings
        set(h_start.ui_popup.lat,'String',column_headers)
        set(h_start.ui_popup.lon,'String',column_headers)
        set(h_start.ui_popup.alt,'String',column_headers)
        set(h_start.ui_popup.yaw,'String',column_headers)
        set(h_start.ui_popup.roll,'String',column_headers)
        set(h_start.ui_popup.pitch,'String',column_headers)
        set(h_start.ui_popup.time,'String',column_headers)
        set(h_start.ui_popup.run,'String',column_headers)
    end
    function ui_util__disable_edit_boxes()
        fields = fieldnames(h_start.ui_edit);
        for k = 1:length(fields)
            set(h_start.ui_edit.(fields{k}),'Enable','off');
        end
        % don't disable ui_text.infile
        set(h_start.ui_edit.infile,'Enable','on');
    end
    function ui_util__enable_edit_boxes()
        fields = fieldnames(h_start.ui_edit);
        for k = 1:length(fields)
            set(h_start.ui_edit.(fields{k}),'Enable','on');
        end
    end
    function ui_util__disable_pushbuttons()
        fields = fieldnames(h_start.ui_dependent_buttons);
        for k = 1:length(fields)
            set(h_start.ui_dependent_buttons.(fields{k}),'Enable','off');
        end
    end
    function ui_util__enable_pushbuttons()
        fields = fieldnames(h_start.ui_dependent_buttons);
        for k = 1:length(fields)
            set(h_start.ui_dependent_buttons.(fields{k}),'Enable','on');
        end
    end

    function ui_util__disable_pushbuttons__write()
        fields = fieldnames(h_start.ui_dependent_buttons);
        for k = 1:length(fields)
            if strfind(lower(fields{k}),'write') == 1
                set(h_start.ui_dependent_buttons.(fields{k}), 'Enable','off');
            end
        end
    end
    function ui_util__enable_pushbuttons__write()
        fields = fieldnames(h_start.ui_dependent_buttons);
        for k = 1:length(fields)
            if strfind(lower(fields{k}),'write') == 1
                set(h_start.ui_dependent_buttons.(fields{k}), 'Enable','on');
            end
        end
    end

    function ui_util__disable_pushbuttons__analyze()
        fields = fieldnames(h_start.ui_dependent_buttons);
        for k = 1:length(fields)
            if strfind(lower(fields{k}),'analyze') == 1
                set(h_start.ui_dependent_buttons.(fields{k}), 'Enable','off');
            end
        end
    end
    function ui_util__enable_pushbuttons__analyze()
        fields = fieldnames(h_start.ui_dependent_buttons);
        for k = 1:length(fields)
            if ~isempty(strfind(lower(fields{k}),'analyze'))
                set(h_start.ui_dependent_buttons.(fields{k}), 'Enable','on');
            end
        end
    end

% (TODO: what happens after 9 platforms?)
    function varargout = fix_TSPI(tspi_class_obj,varargin)
        %   getLLA(this, platform, run_number)
        %   getLLA(this, varargin):
        %       'String-Value' pairs:
        %          'platform'  platform_number (number)
        %          'run'       run_number (number)
        run_num  = [];
        platform = [];
        if(length(varargin) > 0)
            if ischar(varargin{1}) & isempty(str2num(varargin{1}))
                for k = 1:2:length(varargin)
                    switch lower(varargin{k})
                        case 'platform'
                            platform = varargin{k+1};
                        case 'run'
                            run_num = varargin{k+1};
                        case 'preprocess'
                            h_fixtool.platform = platform;
                            h_fixtool.run_num  = run_num;
                            h_fixtool.tspi_class = tspi_class_obj;
                            varargout{1} = fix_TSPI__preprocess();
                            return
                    end
                end
            else
                if length(varargin) == 2
                    platform = varargin{1};
                    run_num  = varargin{2};
                else
                    error('wrong number of inputs')
                end
            end
        else
            platform = 1;
            run_num  = 1;
        end
        
        h_fixtool.platform = platform;
        h_fixtool.run_num  = run_num;
        h_fixtool.tspi_class = tspi_class_obj;
        
        fix_TSPI_util__get_selected_data(platform,run_num);
        
        fix_TSPI__draw_ui()
        
        % whenever the 'String' value in a popup (dropdown) menu is
        % empty, MATLAB issues a warning. We don't care. Does not
        % affect the functionality of the tool. Therefore suppress this
        % warning.
        %    got this string by executing "warning('query','last')"
        %    after I saw the message in the command window.
        warning('off','MATLAB:hg:uicontrol:StringMustBeNonEmpty')
        
        % ===========================================
        % END OF fix_TSPI()
        % ===========================================
        
        function numProblemRegions = fix_TSPI__preprocess()
            numProblemRegions = fix_TSPI__fast_count_problem_regions();
        end
        
        function outstr = fix_TSPI_util__make_figure_title_string()
            outstr = [ ...
                'Platform ' num2str(h_fixtool.platform) ' of ' num2str(tspi_class_obj.nPlats) ...
                ' / Run ' num2str(h_fixtool.run_num) ' of ' num2str(tspi_class_obj.nRuns)];
        end
        
        function fix_TSPI__draw_ui(varargin)
            
            h_fixtool.platName = fix_TSPI_util__make_figure_title_string();
            h_fixtool.hfig_fix = figure('Name',h_fixtool.platName,'DeleteFcn',@fix_TSPI_close_figure,'Tag','tspi_tool');
            set(h_fixtool.hfig_fix,'ButtonDownFcn',@fix_TSPI_cb_bring_into_focus)
            
            screenSize = get(0,'ScreenSize');
            screenSize = screenSize(3:4);
            fixguiSize = [1000 700];
            
            set(h_fixtool.hfig_fix,'Position',[200 screenSize(2)-fixguiSize(2)-100 fixguiSize])
            
            if(h_fixtool.yprFound == 1)
                h_fixtool.hsp(1) = subplot(3,2,1,'Tag','LAT');
                h_fixtool.hsp(2) = subplot(3,2,3,'Tag','LON');
                h_fixtool.hsp(3) = subplot(3,2,5,'Tag','ALT');
                h_fixtool.hsp(4) = subplot(3,2,2,'Tag','YAW'); % YAW
                h_fixtool.hsp(5) = subplot(3,2,4,'Tag','PITCH'); % PITCH
                h_fixtool.hsp(6) = subplot(3,2,6,'Tag','ROLL'); % ROLL
                h_fixtool.select_axes_strings = {'Lat','Lon','Alt','Yaw','Pitch','Roll'};
            else
                h_fixtool.hsp(1) = subplot(3,1,1,'Tag','LAT');
                h_fixtool.hsp(2) = subplot(3,1,2,'Tag','LON');
                h_fixtool.hsp(3) = subplot(3,1,3,'Tag','ALT');
                h_fixtool.select_axes_strings = {'Lat','Lon','Alt'};
            end
            
            linkaxes(h_fixtool.hsp,'x');
            
            h_fixtool.hp_LLA(1) = ...
                plot(h_fixtool.hsp(1), ...
                h_fixtool.latdat, '.-','Tag','LAT');
            
            h_fixtool.hp_LLA(2) = ...
                plot(h_fixtool.hsp(2), ...
                h_fixtool.londat, '.-','Tag','LON');
            
            h_fixtool.hp_LLA(3) = ...
                plot(h_fixtool.hsp(3), ...
                h_fixtool.altdat, '.-','Tag','ALT');
            
            if(h_fixtool.yprFound == 1)
                h_fixtool.hp_YPR(1) = ...
                    plot(h_fixtool.hsp(4), ...
                    h_fixtool.yawdat, '.-','Tag','YAW');
                
                h_fixtool.hp_YPR(2) = ...
                    plot(h_fixtool.hsp(5), ...
                    h_fixtool.pitchdat, '.-','Tag','PITCH');
                
                h_fixtool.hp_YPR(3) = ...
                    plot(h_fixtool.hsp(6), ...
                    h_fixtool.rolldat, '.-','Tag','ROLL');
            end
            
            for kk = 1:length(h_fixtool.hsp)
                hold(h_fixtool.hsp(kk),'on')
            end
            
            
            fix_TSPI_util__identify_problem_regions;
            fix_TSPI_gui_data_status(varargin);
            fix_TSPI_util__update_num_regions;
            % update gui table
            fix_TSPI_util__update_table
            
            xlabel(h_fixtool.hsp(1),'Sample #');ylabel(h_fixtool.hsp(1),'WGS84 Lat(\circ)');%grid on;
            xlabel(h_fixtool.hsp(2),'Sample #');ylabel(h_fixtool.hsp(2),'WGS84 Lon(\circ)');%grid on;
            xlabel(h_fixtool.hsp(3),'Sample #');ylabel(h_fixtool.hsp(3),'WGS84 Alt(m)');%grid on;
            
            if(h_fixtool.yprFound == 1)
                xlabel(h_fixtool.hsp(4),'Sample #');ylabel(h_fixtool.hsp(4),'Yaw (\circ)');%grid on;
                xlabel(h_fixtool.hsp(5),'Sample #');ylabel(h_fixtool.hsp(5),'Pitch (\circ)');%grid on;
                xlabel(h_fixtool.hsp(6),'Sample #');ylabel(h_fixtool.hsp(6),'Roll (\circ)');%grid on;
            end
            
            % Create axes and save handle
            hcmenu_LLA = uicontextmenu(h_fixtool.hfig_fix);
            % Define the context menu items and install their callbacks
            item1 = uimenu(hcmenu_LLA, 'Label', 'Fix Data (this line)',   'Callback',  @fix_TSPI_FIX_IT);
            item3 = uimenu(hcmenu_LLA, 'Label', 'Fix Data (LLA)',         'Callback', {@fix_TSPI_FIX_IT,'LLA'});
            item2 = uimenu(hcmenu_LLA, 'Label', 'Fix Data (LLA and YPR)', 'Callback', {@fix_TSPI_FIX_IT,'all'});
            % Attach the context menu to each line object
            set(h_fixtool.hp_LLA,'uicontextmenu',hcmenu_LLA)
            if(h_fixtool.yprFound == 1)
                hcmenu_YPR = uicontextmenu(h_fixtool.hfig_fix);
                item1 = uimenu(hcmenu_YPR, 'Label', 'Fix Data (this line)',   'Callback',  @fix_TSPI_FIX_IT);
                item3 = uimenu(hcmenu_YPR, 'Label', 'Fix Data (YPR)',         'Callback', {@fix_TSPI_FIX_IT,'YPR'});
                item2 = uimenu(hcmenu_YPR, 'Label', 'Fix Data (LLA and YPR)', 'Callback', {@fix_TSPI_FIX_IT,'all'});
                % Attach the context menu to each lin
                set(h_fixtool.hp_YPR,'uicontextmenu',hcmenu_YPR)
            end
        end
        
        function fix_TSPI_gui_data_status(varargin)
            figpos = get(gcf,'Position');
            
            locatex = figpos(1) + figpos(3)/2;
            locatey = figpos(2) + figpos(4)/2;
            locatex = figpos(1) + figpos(3);
            locatey = figpos(2) + figpos(4);
            
            size_statfig = [250 600];
            
            statfigpos = [locatex locatey-size_statfig(2) size_statfig];
            h_fixtool.hfig_stat = figure('Position',statfigpos, ...
                'ToolBar','none','MenuBar','none', ...
                'DockControls','off', 'NumberTitle', 'off', ...
                'Name','','Tag','tspi_tool');
            
            set(h_fixtool.hfig_stat,'ButtonDownFcn',@fix_TSPI_cb_bring_into_focus)
            
            set(h_fixtool.hfig_stat,'DeleteFcn',@fix_TSPI_close_figure)
            
            figsize = size_statfig;
            
            textheight = 15;
            LRside = 20;
            TopBottom = 20;
            
            for jj = 1:tspi_class_obj.nPlats
                h_fixtool.platStrings{jj} = ['Platform ',num2str(jj)];
            end
            for jj = 1:tspi_class_obj.nRuns
                h_fixtool.runStrings{jj} = ['Run ',num2str(jj)];
            end
            
            kk = 2;
            runplat_space = 15;
            runplat_width = figsize(1)/2 - 2*runplat_space;
            h_fixtool.dropdown_platSelect = uicontrol( ...
                'Style','popupmenu','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight ...
                runplat_width textheight],...
                'String',h_fixtool.platStrings,'Callback',@fix_TSPI_cb_select_platform_and_run, ...
                'Value',h_fixtool.platform);
            
            h_fixtool.dropdown_runSelect = uicontrol( ...
                'Style','popupmenu','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside+runplat_width+runplat_space ...
                figsize(2)-kk*textheight ...
                runplat_width textheight],...
                'String',h_fixtool.runStrings,'Callback',@fix_TSPI_cb_select_platform_and_run, ...
                'Value',h_fixtool.run_num);
            
            % whenever the 'String' value in a popup (dropdown) menu is
            % empty, MATLAB issues a warning. We don't care. Does not
            % affect the functionality of the tool. Therefore suppress this
            % warning.
            %    got this string by executing "warning('query','last')"
            %    after I saw the message in the command window.
            warning('off','MATLAB:hg:uicontrol:StringMustBeNonEmpty')
            
            kk = kk + 2.0;
            h_fixtool.dropdown_axesSelect = uicontrol( ...
                'Style','popupmenu','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'String',h_fixtool.select_axes_strings,'Callback',@fix_TSPI_cb_select_axes);
            h_fixtool.select_axes = h_fixtool.hsp(1);
            
            kk = kk + 2.0;
            %                 select_regions_strings = num2cell(get_num_regions_in_axes(h_fixtool.select_axes));
            select_regions_list = 1:fix_TSPI_util__get_num_regions_in_axes(h_fixtool.select_axes);
            h_fixtool.dropdown_regionsSelect = uicontrol( ...
                'Style','popupmenu','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'String',select_regions_list,'Callback',@fix_TSPI_cb_select_region);
            
            % run the callbacks for the select_axes and select_region
            % dropdowns
            fix_TSPI_cb_select_axes(h_fixtool.dropdown_axesSelect);
            fix_TSPI_cb_select_region(h_fixtool.dropdown_regionsSelect);
            
            kk = kk + 2.5;
            prevnext_width = 40;
            goto_width     = 100;
            button_spacing = 15;
            h_fixtool.button_goto_region = uicontrol( ...
                'Style','pushbutton','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside+prevnext_width+1*button_spacing ...
                figsize(2)-kk*textheight goto_width textheight+10],...
                'String','Go to region', ...
                'Callback',@fix_TSPI_cb_goto_region);
            
            h_fixtool.button_goto_region_next = uicontrol( ...
                'Style','pushbutton','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside+prevnext_width+goto_width+2*button_spacing ...
                figsize(2)-kk*textheight prevnext_width textheight+10],...
                'String','Next', ...
                'Callback',{@fix_TSPI_cb_goto_region_next,'direction',1});
            h_fixtool.button_goto_region_prev = uicontrol( ...
                'Style','pushbutton','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight prevnext_width textheight+10],...
                'String','Prev', ...
                'Callback',{@fix_TSPI_cb_goto_region_next,'direction',-1});
            
            kk = kk + 2.0;
            h_fixtool.button_restore_full_view = uicontrol( ...
                'Style','pushbutton','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight size_statfig(1)-2*LRside textheight+10],...
                'String','Restore full view', ...
                'Callback',@fix_TSPI_cb_restore_full_view);
            %                 [LRside 2.0*TopBottom 120 textheight]
            
            kk = kk + 1.5;
            h_fixtool.num_regions_lat = uicontrol( ...
                'Style','text','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'HorizontalAlignment','left', ...
                'String',['Regions in Lat: ',num2str(fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(1)))]);
            kk = kk + 1;
            h_fixtool.num_regions_lon = uicontrol( ...
                'Style','text','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'HorizontalAlignment','left', ...
                'String',['Regions in Lon: ',num2str(fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(2)))]);
            kk = kk + 1;
            h_fixtool.num_regions_alt = uicontrol( ...
                'Style','text','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'HorizontalAlignment','left', ...
                'String',['Regions in Alt: ',num2str(fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(3)))]);
            kk = kk + 1;
            h_fixtool.num_regions_yaw = uicontrol( ...
                'Style','text','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'HorizontalAlignment','left', ...
                'String',['Regions in Yaw: ']);
            kk = kk + 1;
            h_fixtool.num_regions_pit = uicontrol( ...
                'Style','text','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'HorizontalAlignment','left', ...
                'String',['Regions in Pitch: ']);
            kk = kk + 1;
            h_fixtool.num_regions_rol = uicontrol( ...
                'Style','text','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',[LRside figsize(2)-kk*textheight figsize(1)-2*LRside textheight],...
                'HorizontalAlignment','left', ...
                'String',['Regions in Roll: ']);
            
            h_fixtool.xlim_full_view = xlim(h_fixtool.hsp(1));

            kk = kk + 4;
            h_fixtool.ui_table = uitable('Parent',h_fixtool.hfig_stat, ...
                'Position', [LRside TopBottom/2+textheight+30 size_statfig(1)-2*LRside size_statfig(2)-kk*textheight], ...
                'Data', []);
            
            for krun = 1:h_fixtool.tspi_class.nRuns
                table_rownames{krun} = ['Run ',num2str(krun)];
            end
            for kplat = 1:h_fixtool.tspi_class.nPlats
                table_colnames{kplat} = ['Plat ',num2str(kplat)];
            end
            set(h_fixtool.ui_table,'RowName',table_rownames);
            set(h_fixtool.ui_table,'ColumnName',table_colnames);
            set(h_fixtool.ui_table,'ColumnWidth',num2cell(50*ones(1,h_fixtool.tspi_class.nRuns)));
            set(h_fixtool.ui_table,'TooltipString','click to go to data')
            set(h_fixtool.ui_table,'CellSelectionCallback',@fix_TSPI_cb_goto_selected_table_point)
            fix_TSPI_util__update_table;
            
            button_finished_size = 200;
            %             button_finished_pos = [LRside TopBottom/2 80 textheight+10];
            button_finished_pos = [size_statfig(1)/2-button_finished_size/2 TopBottom/2 button_finished_size textheight+10];
            h_fixtool.button_write_file = uicontrol( ...
                'Style','pushbutton','Parent',h_fixtool.hfig_stat, ...
                'Units','Pixels', ...
                'Position',button_finished_pos,...
                'String','Finished!', 'FontSize', 12, ...
                'Callback',@fix_TSPI__finalize_and_return_to_main_tool);
        end
        
        
        function fix_TSPI_cb_goto_selected_table_point(varargin)
            tableData            = varargin{1};
            cellSelectChangeData = varargin{2};
            
            selectedInd = cellSelectChangeData.Indices;
            
            % get only the first select table point (if multiples)
            if ~isempty(selectedInd)
                kPlat = selectedInd(1,2);
                kRun  = selectedInd(1,1);
            else
                return
            end
            
            % update gui dropdowns, then get this new data
            fix_TSPI_util__set_dropdown_platform_and_run(kPlat,kRun)
            fix_TSPI_cb_select_platform_and_run;
            
        end
        
        function fix_TSPI_util__set_dropdown_platform_and_run(kPlat,kRun)
            set(h_fixtool.dropdown_platSelect,'Value',kPlat)
            set(h_fixtool.dropdown_runSelect, 'Value',kRun)
        end
        
        function fix_TSPI_util__update_table(varargin)
            fix_TSPI__fast_count_problem_regions()
            set(h_fixtool.ui_table,'Data',h_fixtool.numTestsFailed)
            if(all(h_fixtool.numTestsFailed == 0))
                ui_util__enable_pushbuttons__write();
            else
                ui_util__disable_pushbuttons__write();
            end
        end
        
        function fix_TSPI_cb_select_platform_and_run(varargin)
            kPlat = get(h_fixtool.dropdown_platSelect,'Value');
            kRun  = get(h_fixtool.dropdown_runSelect, 'Value');
            
            fix_TSPI_util__get_selected_data(kPlat,kRun);
            
            set(h_fixtool.hfig_fix,'Name',fix_TSPI_util__make_figure_title_string());
            
            fix_TSPI_util__update_plotted_data();
            
            % now clear all existing patches
            for jj = 1:length(h_fixtool.hsp)
                hpatches = findobj(get(h_fixtool.hsp(jj),'Children'),'Type','patch');
                delete(hpatches)
            end
            
            fix_TSPI_util__identify_problem_regions;
            fix_TSPI_util__update_num_regions;
            % update gui table
            fix_TSPI_util__update_table
            
        end
        function fix_TSPI_cb_restore_full_view(varargin)
            set(h_fixtool.hsp(1),'xlim',h_fixtool.xlim_full_view)
        end
        function fix_TSPI_cb_goto_region(varargin)
            
            % make sure that the latest is grabbed
            fix_TSPI_cb_select_axes(h_fixtool.dropdown_axesSelect)
            fix_TSPI_cb_select_region(h_fixtool.dropdown_regionsSelect)
            
            % get region width
            
            xdat = get(h_fixtool.select_region,'xdata');
            wx = max(xdat) - min(xdat);
            wx = max(wx,50);
            cx = 0.5*(min(xdat) + max(xdat));
            
            xlim(h_fixtool.select_axes, cx + 1*[-wx wx]);
            
        end
        function fix_TSPI_cb_goto_region_next(varargin)
            
            for k = 1:2:nargin
                switch lower(varargin{k})
                    case 'direction'
                        increment = varargin{k+1};
                end
            end
            
            % increment region number (if region is last in this axes,
            % then increment to next axes)
            maxVal_axes   = length(get(h_fixtool.dropdown_axesSelect,'String'));
            maxVal_region = length(get(h_fixtool.dropdown_regionsSelect,'String'));
            
            kNextRegion = get(h_fixtool.dropdown_regionsSelect,'Value') + increment;
            kNextAxes   = get(h_fixtool.dropdown_axesSelect,'Value') + increment;
            
            if kNextAxes > maxVal_axes
                kNextAxes = 1;
            elseif kNextAxes == 0
                kNextAxes = maxVal_axes;
            end
            
            hNextAxes = h_fixtool.hsp(kNextAxes);
            nRegionsInNextAxes = fix_TSPI_util__get_num_regions_in_axes(hNextAxes);
            
            if kNextRegion > maxVal_region
                % WRAPPING AROUND FROM LAST AXES
                set(h_fixtool.dropdown_axesSelect,'Value',kNextAxes)
                set(h_fixtool.dropdown_regionsSelect,'Value',1)
            elseif kNextRegion == 0
                % WRAPPING AROUND FROM FIRST AXES
                set(h_fixtool.dropdown_axesSelect,'Value',kNextAxes)
                set(h_fixtool.dropdown_regionsSelect,'Value',nRegionsInNextAxes)
            else
                set(h_fixtool.dropdown_regionsSelect,'Value',kNextRegion)
            end
            
            % make sure that the latest is grabbed
            fix_TSPI_cb_select_axes(h_fixtool.dropdown_axesSelect)
            fix_TSPI_cb_select_region(h_fixtool.dropdown_regionsSelect)
            
            % get region width
            if ~isempty(h_fixtool.select_region)
                xdat = get(h_fixtool.select_region,'xdata');
                wx = max(xdat) - min(xdat);
                wx = max(wx,50);
                cx = 0.5*(min(xdat) + max(xdat));
                
                xlim(h_fixtool.select_axes, cx + 1*[-wx wx]);
            end
        end
        function fix_TSPI_cb_select_axes(varargin)
            hdrop = varargin{1};
            h_fixtool.select_axes = h_fixtool.hsp(get(hdrop,'Value'));
            nRegions = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.select_axes);
            set(h_fixtool.dropdown_regionsSelect,'String',1:nRegions)
        end
        function fix_TSPI_cb_select_region(varargin)
            hdrop = varargin{1};
            hpatches = fix_TSPI_util__get_patches_in_axes(h_fixtool.select_axes);
            if ~isempty(hpatches)
                h_fixtool.select_region = hpatches(get(hdrop,'Value'));
            else
                % if no problem regions, set to empty
                h_fixtool.select_region = [];
            end
        end
        
        
        
        
        function varargout = fix_TSPI_util__identify_problem_regions()
            
            plot_flag = 1;
            nAxes = length(h_fixtool.hsp);
            numProblemRegions = 0;
            for isp = 1:length(h_fixtool.hsp)
                switch isp
                    case 1
                        dat = h_fixtool.latdat;
                        hax = h_fixtool.hsp(1);
                    case 2
                        dat = h_fixtool.londat;
                        hax = h_fixtool.hsp(2);
                    case 3
                        dat = h_fixtool.altdat;
                        hax = h_fixtool.hsp(3);
                    case 4
                        dat = h_fixtool.yawdat;
                        hax = h_fixtool.hsp(4);
                    case 5
                        dat = h_fixtool.pitchdat;
                        hax = h_fixtool.hsp(5);
                    case 6
                        dat = h_fixtool.rolldat;
                        hax = h_fixtool.hsp(6);
                end
                
                % get any previous patches and DELETE THEM
                hpatches = findobj(get(hax,'Children'),'Type','patch');
                delete(hpatches)
                
                test_1_pass = fix_TSPI_region_test_1(hax,dat,plot_flag);
                test_2_pass = fix_TSPI_region_test_2(hax,    plot_flag);
                
                if(nargout > 0)
                    numProblemRegions = numProblemRegions + fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(isp));
                end
            end
            
            if(nargout > 0)
                varargout{1} = numProblemRegions;
            end
            
        end
        function num_regions = fix_TSPI_util__get_num_regions_in_axes(hax)
            num_regions = length(findobj(get(hax,'Children'),'Type','patch'));
        end
        function fix_TSPI_util__update_plotted_data()
            
            % delete all "selected_points" line plots
            htmp1 = findobj(get(h_fixtool.hsp(1),'Children'),'Tag','selected_points');
            htmp2 = findobj(get(h_fixtool.hsp(2),'Children'),'Tag','selected_points');
            htmp3 = findobj(get(h_fixtool.hsp(3),'Children'),'Tag','selected_points');
            delete(htmp1)
            delete(htmp2)
            delete(htmp3)
            
            set(h_fixtool.hp_LLA(1), ...
                'xdata',1:length(h_fixtool.latdat), ...
                'ydata',h_fixtool.latdat);
            set(h_fixtool.hp_LLA(2), ...
                'xdata',1:length(h_fixtool.londat), ...
                'ydata',h_fixtool.londat);
            set(h_fixtool.hp_LLA(3), ...
                'xdata',1:length(h_fixtool.altdat), ...
                'ydata',h_fixtool.altdat);
            
            nSamples = length(h_fixtool.latdat);
            xlim(h_fixtool.hsp(1),[1 nSamples])
            xlim(h_fixtool.hsp(2),[1 nSamples])
            xlim(h_fixtool.hsp(3),[1 nSamples])
            
            if(h_fixtool.yprFound == 1)
                % delete all "selected_points" line plots
                htmp1 = findobj(get(h_fixtool.hsp(4),'Children'),'Tag','selected_points');
                htmp2 = findobj(get(h_fixtool.hsp(5),'Children'),'Tag','selected_points');
                htmp3 = findobj(get(h_fixtool.hsp(6),'Children'),'Tag','selected_points');
                delete(htmp1)
                delete(htmp2)
                delete(htmp3)
                
                set(h_fixtool.hsp(4:6),'Visible','on')
                hchild = get(h_fixtool.hsp(4),'Children');
                set(hchild,'Visible','on');
                hchild = get(h_fixtool.hsp(5),'Children');
                set(hchild,'Visible','on');
                hchild = get(h_fixtool.hsp(6),'Children');
                set(hchild,'Visible','on');
                set(h_fixtool.hp_YPR(1), ...
                    'xdata',1:length(h_fixtool.yawdat), ...
                    'ydata',h_fixtool.yawdat);
                set(h_fixtool.hp_YPR(2), ...
                    'xdata',1:length(h_fixtool.pitchdat), ...
                    'ydata',h_fixtool.pitchdat);
                set(h_fixtool.hp_YPR(3), ...
                    'xdata',1:length(h_fixtool.rolldat), ...
                    'ydata',h_fixtool.rolldat);
                nSamples = length(h_fixtool.yawdat);
                xlim(h_fixtool.hsp(4),[1 nSamples])
                xlim(h_fixtool.hsp(5),[1 nSamples])
                xlim(h_fixtool.hsp(6),[1 nSamples])
            else
                % if there are axes for previous YPR data, then clear
                % the data in the handles
                if(isfield(h_fixtool,'hp_YPR'))
                    hpatches = findobj(get(h_fixtool.hsp(4),'Children'),'Type','patch');
                    delete(hpatches);
                    hpatches = findobj(get(h_fixtool.hsp(5),'Children'),'Type','patch');
                    delete(hpatches);
                    hpatches = findobj(get(h_fixtool.hsp(6),'Children'),'Type','patch');
                    delete(hpatches);
                    set(h_fixtool.hsp(4:6),'Visible','off')
                    set(h_fixtool.hp_YPR(1), ...
                        'xdata',[], ...
                        'ydata',[]);
                    set(h_fixtool.hp_YPR(2), ...
                        'xdata',[], ...
                        'ydata',[]);
                    set(h_fixtool.hp_YPR(3), ...
                        'xdata',[], ...
                        'ydata',[]);
                end
            end
            
        end
        function fix_TSPI_util__get_selected_data(platform,run_num)
            h_fixtool.platform = platform;
            h_fixtool.run_num  = run_num;
            
            LLADAT = tspi_class_obj.getLLA('platform',platform,'run',run_num);
            h_fixtool.latdat = LLADAT(:,1);
            h_fixtool.londat = LLADAT(:,2);
            h_fixtool.altdat = LLADAT(:,3);
            
            r2d = 180/pi;
            try
                YPRDAT = r2d * tspi_class_obj.getYPR('platform',platform,'run',run_num);
                h_fixtool.yawdat   = YPRDAT(:,1);
                h_fixtool.pitchdat = YPRDAT(:,2);
                h_fixtool.rolldat   = YPRDAT(:,3);
                h_fixtool.yprFound = 1;
            catch
                if(isfield(h_fixtool,'yawdat'))
                    rmfield(h_fixtool,'yawdat')
                    rmfield(h_fixtool,'pitchdat')
                    rmfield(h_fixtool,'rolldat')
                end
                h_fixtool.yprFound = 0;
            end
        end
        function hpatches = fix_TSPI_util__get_patches_in_axes(hax)
            hpatches = findobj(get(hax,'Children'),'Type','patch');
            % flip direction (last patch added is the first one in the
            % list)
            hpatches = hpatches(end:-1:1);
        end
        function fix_TSPI_util__update_num_regions(varargin)
            % compute the remaining regions not addressed
            nregions(1) = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(1));
            nregions(2) = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(2));
            nregions(3) = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(3));
            
            if(nregions(1) == 0); tag_str{1} = ' FINISHED!'; else; tag_str{1} = ''; end
            if(nregions(2) == 0); tag_str{2} = ' FINISHED!'; else; tag_str{2} = ''; end
            if(nregions(3) == 0); tag_str{3} = ' FINISHED!'; else; tag_str{3} = ''; end
            
            set(h_fixtool.num_regions_lat,'String',['Regions in Lat: ' num2str(nregions(1)) tag_str{1}])
            set(h_fixtool.num_regions_lon,'String',['Regions in Lon: ' num2str(nregions(2)) tag_str{2}])
            set(h_fixtool.num_regions_alt,'String',['Regions in Alt: ' num2str(nregions(3)) tag_str{3}])
            
            if(h_fixtool.yprFound == 1)
                nregions(4) = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(4));
                nregions(5) = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(5));
                nregions(6) = fix_TSPI_util__get_num_regions_in_axes(h_fixtool.hsp(6));
                if(nregions(4) == 0); tag_str{4} = ' FINISHED!'; else; tag_str{4} = ''; end
                if(nregions(5) == 0); tag_str{5} = ' FINISHED!'; else; tag_str{5} = ''; end
                if(nregions(6) == 0); tag_str{6} = ' FINISHED!'; else; tag_str{6} = ''; end
                set(h_fixtool.num_regions_yaw,'String',['Regions in Yaw: '   num2str(nregions(4)) tag_str{4}])
                set(h_fixtool.num_regions_pit,'String',['Regions in Pitch: ' num2str(nregions(5)) tag_str{5}])
                set(h_fixtool.num_regions_rol,'String',['Regions in Roll: '  num2str(nregions(6)) tag_str{6}])
            else
                set(h_fixtool.num_regions_yaw,'String','')
                set(h_fixtool.num_regions_pit,'String','')
                set(h_fixtool.num_regions_rol,'String','')
            end
            
            set(h_fixtool.dropdown_regionsSelect,'String',1:fix_TSPI_util__get_num_regions_in_axes(h_fixtool.select_axes))
            
        end
        
        function [center,width,startstop,iregion] = label_regions(logical_test)
            
            iregion = zeros(size(logical_test));
            
            iregcnt = 0;
            % if first is true, then start "new region"
            if(logical_test(1))
                iregcnt = iregcnt + 1;
                iregion(1) = 1;
            end
            
            for k = 2:length(logical_test)
                % is detection?
                if(logical_test(k))
                    % is new region? look to previous value
                    % part of existing region
                    if(iregion(k-1) ~= 0)
                        iregion(k) = iregion(k-1);
                    else
                        iregcnt = iregcnt + 1;
                        iregion(k) = iregcnt;
                    end
                end
            end
            
            numRegions = iregcnt;
            
            center = zeros(1,numRegions);
            width  = center;
            startstop = zeros(numRegions,2);
            
            idx = 1:length(logical_test);
            for region_label = 1:numRegions
                %
                isolated_region     = (iregion == region_label);
                idx_isolated_region = idx(isolated_region);
                center(region_label) = mean(idx_isolated_region);
                width(region_label) = sum(isolated_region*1);
                startstop(region_label,1) = find(isolated_region,1,'first');
                startstop(region_label,2) = find(isolated_region,1,'last');
            end
            
        end
        
        function varargout = fix_TSPI__fast_count_problem_regions()
            % using selected data in h_fixtool.latdat, londat, altdat,
            % etc..., run the tests WITHOUT plotting
            numTestsFailed = zeros(h_fixtool.tspi_class.nRuns, h_fixtool.tspi_class.nPlats);
            % save off the plat/run numbers currently selected since
            curPlat = h_fixtool.platform;
            curRun  = h_fixtool.run_num;
            for kplat = 1:h_fixtool.tspi_class.nPlats
                for krun = 1:h_fixtool.tspi_class.nRuns
                    fix_TSPI_util__get_selected_data(kplat,krun);
                    test_1_res_lat   = fix_TSPI_region_test_1([],'lat',0);
                    test_1_res_lon   = fix_TSPI_region_test_1([],'lon',0);
                    test_1_res_alt   = fix_TSPI_region_test_1([],'alt',0);
                    test_1_res_yaw   = fix_TSPI_region_test_1([],'yaw',0);
                    test_1_res_pitch = fix_TSPI_region_test_1([],'pitch',0);
                    test_1_res_roll  = fix_TSPI_region_test_1([],'roll',0);
                    test_2_res_lla   = fix_TSPI_region_test_2('lla',0);
                    test_2_res_yrp   = fix_TSPI_region_test_2('ypr',0);
                    numTestsFailed(krun,kplat) = sum(1.0*(~[ ...
                        test_1_res_lat test_1_res_lon test_1_res_alt ...
                        test_1_res_yaw test_1_res_pitch test_1_res_roll ...
                        test_2_res_lla test_2_res_yrp]));
                end
            end
            % reset the platform/run numbers
            h_fixtool.platform = curPlat;
            h_fixtool.run_num  = curRun;
            
            h_fixtool.numTestsFailed = numTestsFailed;
            
            if(nargout > 0)
                varargout{1} = numTestsFailed;
            end
        end
        
        function test_passed = fix_TSPI_region_test_1(hax,dat,plot_flag,varargin)
            % set plot_flag = 1  if first pass and want to draw the
            % patches,
            % set plot_flag = 0  if only want to perform the tests.
            if(nargin >= 4)
                idx = varargin{1};
            else
                idx = 1:length(h_fixtool.latdat);
            end
            
            if(isempty(dat))
                return
            end
            
            if ischar(dat)
                switch lower(dat)
                    case 'lat'
                        dat = h_fixtool.latdat;
                    case 'lon'
                        dat = h_fixtool.londat;
                    case 'alt'
                        dat = h_fixtool.altdat;
                    case 'yaw'
                        dat = h_fixtool.yawdat;
                    case 'pitch'
                        dat = h_fixtool.pitchdat;
                    case 'roll'
                        dat = h_fixtool.rolldat;
                    otherwise
                        error(['unrecognized data "', dat, '"'])
                end
            end
            testdat = isnan(dat(idx));
            
            if(length(testdat) > 0)
                [center,width,startstop,iregion] = label_regions(testdat);
            else
                test_passed = logical(1);
                return
            end
            % DRAW PATCHES OVER REGIONS OF INTEREST
            nRegions = length(center);
            test_passed = logical(1);
            nDataPoints = length(testdat);
            
            if(nRegions > 0)
                if(nargin == 5)
                    regionsToTest = varargin{2};
                else
                    regionsToTest = 1:nRegions;
                end
                test_passed = logical(zeros(1,nRegions));
                for ireg = regionsToTest
                    if(plot_flag == 0)
                        continue
                    end
                    
                    % if the end of this region is the last sample
                    % point then DON'T make this a problem area.
                    if(startstop(ireg,2) == nDataPoints)
                        continue
                    end
                    ylims = ylim(hax);
                    px = [startstop(ireg,1)-0.5 startstop(ireg,2)+0.5 startstop(ireg,2)+0.5 startstop(ireg,1)-0.5];
                    py = [ylims(1) ylims(1) ylims(2) ylims(2)];
                    mx = 0.5*(min(px)+max(px));
                    my = 0.5*(min(py)+max(py));
                    hpatch = patch(hax,px,py,'k');
                    set(hpatch,'FaceAlpha',0.3,'EdgeColor','none')
                end
            end
        end
        function test_passed = fix_TSPI_region_test_2(hax,plot_flag,varargin)
            % find regions with POTENTIAL STALE DATA (REPEATED VALUES
            % ON ALL THREE LINES)
            if(nargin >= 3)
                idx = varargin{1};
            else
                idx = 1:length(h_fixtool.latdat);
            end
            
            if(ishandle(hax))
                if any(hax == [h_fixtool.hsp(1) h_fixtool.hsp(2) h_fixtool.hsp(3)])
                    ddat1 = diff(h_fixtool.latdat(idx));
                    ddat2 = diff(h_fixtool.londat(idx));
                    ddat3 = diff(h_fixtool.altdat(idx));
                elseif any(hax == [h_fixtool.hsp(4) h_fixtool.hsp(5) h_fixtool.hsp(6)])
                    ddat1 = diff(h_fixtool.yawdat(idx));
                    ddat2 = diff(h_fixtool.pitchdat(idx));
                    ddat3 = diff(h_fixtool.rolldat(idx));
                end
            elseif(ischar(hax))
                switch lower(hax)
                    case 'lla'
                        ddat1 = diff(h_fixtool.latdat(idx));
                        ddat2 = diff(h_fixtool.londat(idx));
                        ddat3 = diff(h_fixtool.altdat(idx));
                    case 'ypr'
                        ddat1 = diff(h_fixtool.yawdat(idx));
                        ddat2 = diff(h_fixtool.pitchdat(idx));
                        ddat3 = diff(h_fixtool.rolldat(idx));
                end
            end
            
            % consecutive stale data threshold
            h_fixtool.staleDataThreshold = 3;
            
            testdat = (ddat1 == 0) & (ddat2 == 0) & (ddat3 == 0);
            
            if length(testdat) > 0
                [center,width,startstop,iregion] = label_regions(testdat);
                if(length(center) == 0)
                    test_passed = logical([1 1 1]);
                    return
                else
                    % need to perform some corrections due to diff()
                    % operation on the original data.
                    width = width + 1;
                    startstop(:,2) = startstop(:,2) + 1;
                end
            else
                test_passed = logical([1 1 1]);
                return
            end
            
            %                     [center,width,startstop,iregion] = label_regions(testdat);
            % DRAW PATCHES OVER REGIONS OF INTEREST
            nRegions = length(center);
            iproblem = 0;
            
            if(nargin == 4)
                regionsToTest = varargin{2};
            else
                regionsToTest = 1:nRegions;
            end
            
            %             test_passed = logical(ones(1,3*nRegions));
            ireg_failed = 0;
            for ireg = regionsToTest
                if(width(ireg) >= h_fixtool.staleDataThreshold)
                    ireg_failed = ireg_failed + 1;
                    %                     test_passed(ireg) = logical(0);
                    if(plot_flag == 0)
                        continue
                    end
                    iproblem = iproblem + 1;
                    ylims = ylim(hax);
                    px = [startstop(ireg,1)-0.5 startstop(ireg,2)+0.5 startstop(ireg,2)+0.5 startstop(ireg,1)-0.5];
                    py = [ylims(1) ylims(1) ylims(2) ylims(2)];
                    hpatch = patch(hax,px,py,'r');
                    set(hpatch,'FaceAlpha',0.5,'EdgeColor','none')
                    %                             hpatches(iproblem) = hpatch;
                end
            end
            test_passed = logical(zeros(1,3*ireg_failed));
        end
        
        function fix_TSPI_FIX_IT(varargin)
            % since we have to click on the line to access the context menu
            % and click on the item that calls this callback, the current
            % object (gco) is the line itself!
            option_fix_str = 'single';
            for kk = 1:nargin
                switch lower(varargin{kk})
                    case 'all'
                        option_fix_str = varargin{kk};
                    case 'lla'
                        option_fix_str = varargin{kk};
                    case 'ypr'
                        if(h_fixtool.yprFound == 1)
                            option_fix_str = varargin{kk};
                        else
                            error('tried to fix YPR when no YPR is present')
                        end
                end
            end
            curve = gco;
            curvehax = get(curve,'Parent');
            if(length(curve) == 1)
                set(curve,'UserData','master');
            else
                curve = findobj(curve,'UserData','master');
            end
            
            xd=get(curve,'XData');
            yd=get(curve,'YData');
            
            if ~isprop(h_fixtool.hfig_fix,'KeyPressFcn')
                interact_method_used = 1;
                % METHOD #1 (for pre-MATLAB R2016)
                uiwait(msgbox({'Click data points to keep';'Press [ENTER] when done'},'Fix-It'));
                gleft = ginput();
                numClickedPoints = size(gleft,1);
                
                % now re-pack to nearest data
                xpick = gleft(:,1);
                ypick = gleft(:,2);
            else
                interact_method_used = 2;
                % METHOD #2 (a little more dynamic and helpful)
                %   wait for enter key to be pressed
                set(curve,'ButtonDownFcn',@fix_TSPI_cb_click_points_for_spline)
                set(h_fixtool.hfig_fix,'KeyPressFcn',{@fix_TSPI_cb_click_points_for_spline_END,curvehax})
                uiwait
                set(hfig,'KeyPressFcn',[])
                
                hpsel = findobj(get(curvehax,'Children'),'Tag','selected_points');
                numClickedPoints = length(hpsel);
                xpick = zeros(numClickedPoints,1);
                ypick = zeros(numClickedPoints,1);
                for k = 1:numClickedPoints
                    xpick(k,1) = get(hpsel(k),'XData');
                    ypick(k,1) = get(hpsel(k),'YData');
                end
            end
            
            xvals = [];
            xidxs = [];
            yvals = [];
            
            for k = 1 : numClickedPoints
                d = (xpick(k)-xd).^2 + (ypick(k)-yd).^2;
                [m,midx]=min(d);
                % if the point has already been selected, skip
                if(isempty(find(xidxs == midx)))
                    xidxs = [xidxs midx];
                    hold on;
                end
            end
            
            % determine the number of curves to update
            switch lower(option_fix_str)
                case 'single'
                    curves = curve;
                    curvehaxes = curvehax;
                case 'all'
                    curves = [h_fixtool.hp_LLA h_fixtool.hp_YPR];
                    curvehaxes = h_fixtool.hsp;
                case 'lla'
                    curves = h_fixtool.hp_LLA;
                    curvehaxes = h_fixtool.hsp;
                case 'ypr'
                    curves = h_fixtool.hp_YPR;
                    curvehaxes = h_fixtool.hsp;
            end
            
            % perform fixes on all requested data curves
            for kk = 1:length(curves)
                
                thiscurve = curves(kk);
                xd = get(thiscurve,'XData');
                yd = get(thiscurve,'YData');
                
                % sort xidxs
                [xidxs,sortidx] = sort(xidxs);
                if interact_method_used == 1
                    tmphp = plot(curvehaxes(kk),xd(xidxs),yd(xidxs),'ro');
                    set(tmphp,'Tag','selected_points');
                end
                
                % compute new data based on spline interpolation
                pp = spline(xidxs,yd(xidxs), (xidxs(1)+1) : (xidxs(end)-1) );
                yd( (xidxs(1)+1) : (xidxs(end)-1) ) = pp;
                
                set(thiscurve,'YData',yd);
            end
            
            % after fixing a region, update the internal data table
            % with the new values (keeping the table up-to-date for
            % when the "write to file" button is pressed).
            fix_TSPI_util__update_internal_data_table();
            
            % re-analyze data set for problem areas
            fix_TSPI_util__identify_problem_regions
            
            % update mini gui
            fix_TSPI_util__update_num_regions
            
            % update gui table
            fix_TSPI_util__update_table
            
            % delete selected points!
            delete(hpsel)
            
            % remove ButtonDownFcn callback on the data line
            set(curve,'ButtonDownFcn',[])
            
            function fix_TSPI_cb_click_points_for_spline(varargin)
                hp      = varargin{1};
                hit_obj = varargin{2};
                
                % axes parent of plot
                hax = get(hp,'Parent');
                
                % get current axes point
                cpt = get(hax,'CurrentPoint');
                
                % find nearest point in plot
                xdat = get(hp,'XData');
                ydat = get(hp,'YData');
                
                xpt = cpt(1,1);
                ypt = cpt(1,2);
                
                dist = sqrt((xdat - xpt).^2 + (ydat - ypt).^2);
                [mindist,minidx] = min(dist(:));
                
                % plot selected point
                xsel = xdat(minidx);
                ysel = ydat(minidx);
                hold(hax,'on')
                hpsel = plot(xsel,ysel,'ro','ButtonDownFcn',@fix_TSPI_cb_click_on_selected);
                set(hpsel,'Tag','selected_points');
            end
            
            function fix_TSPI_cb_click_on_selected(varargin)
                hpselected = varargin{1};
                delete(hpselected)
            end
            
            function fix_TSPI_cb_click_points_for_spline_END(varargin)
                hfig = varargin{1};
                keydat_obj = varargin{2};
                
                if strcmpi(keydat_obj.Key,'return')
                    uiresume
                end
            end
            
            function fix_TSPI_reanalyze_data_for_problems(hline)
                
                hax = get(hline,'Parent');
                
                xdat = get(hline,'XData');
                ydat = get(hline,'YData');
                
                hpatches = findobj(get(hax,'Children'),'Type','patch');
                
                % delete all patches and redo analysis (which highlights
                % any areas that are still problems)
                delete(hpatches)
                
                fix_TSPI_util__identify_problem_regions
                
            end
        end
        
        function fix_TSPI_util__update_internal_data_table()
            % get the latest data
            h_fixtool.latdat = get(h_fixtool.hp_LLA(1),'YData');
            h_fixtool.londat = get(h_fixtool.hp_LLA(2),'YData');
            h_fixtool.altdat = get(h_fixtool.hp_LLA(3),'YData');
            
            % set latest data to table
            setLLA(h_fixtool.tspi_class,[h_fixtool.latdat(:) h_fixtool.londat(:) h_fixtool.altdat(:)], ...
                'platform',h_fixtool.platform,'run',h_fixtool.run_num);
            
            if(h_fixtool.yprFound == 1)
                h_fixtool.yawdat   = get(h_fixtool.hp_YPR(1),'YData');
                h_fixtool.pitchdat = get(h_fixtool.hp_YPR(2),'YData');
                h_fixtool.rolldat  = get(h_fixtool.hp_YPR(3),'YData');
                setYPR(h_fixtool.tspi_class,[h_fixtool.yawdat(:) h_fixtool.pitchdat(:) h_fixtool.rolldat(:)], ...
                    'platform',h_fixtool.platform,'run',h_fixtool.run_num);
            end
        end
        
        function fix_TSPI__finalize_and_return_to_main_tool(varargin)
            % need to put modified data back into the table
            % "tspi_class.tableData"
            fix_TSPI_util__update_internal_data_table();
            
            fix_TSPI_close_figure;
            
            % make sure to see if all problems have been solved
            fix_TSPI__fast_count_problem_regions();
            if(h_fixtool.numTestsFailed == 0)
                ui_util__enable_pushbuttons__write();
            else
                ui_util__disable_pushbuttons__write();
            end
        end
        
        function fix_TSPI_cb_bring_into_focus(varargin)
            figure(h_fixtool.hfig_stat)
            figure(h_fixtool.hfig_fix)
        end
        
        function fix_TSPI_close_figure(varargin)
            delete(h_fixtool.hfig_stat)
            delete(h_fixtool.hfig_fix)
            warning('on')
        end
    end

    function ui_cb__close_main_gui(varargin)
        warning('on')
        myObjs = findobj('Tag','tspi_tool');
        delete(myObjs);
    end
end
