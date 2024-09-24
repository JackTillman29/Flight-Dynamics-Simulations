function cb_WaveFftStruct(timeWindow)

if ~exist('WaveObj')
    error('needs "WaveObj" (located in MATLAB\Analysis)')
end
if ~exist('rectwin')
    warning('may not have access to window functions (located in MATLAB\Windows)')
end

% get figure handle
hfig = gcf;

% get axes handle
hax = findobj('Parent',hfig,'Type','Axes');

% get data
hchild = get(hax(1),'Children');

% analyze each waveform on axes
for i = 1:length(hchild)
    
    x = get(hchild(i),'XData');
    y = get(hchild(i),'YData');
    
    if(nargin > 0)
        tlim = get(gca,'XLim');
        itStart = find(x >= tlim(1),1);
        itEnd   = find(x >= tlim(2),1);
        if(isempty(itStart))
            itStart = 1;
        end
        if(isempty(itEnd))
            itEnd = length(x);
        end
        x = x(itStart:itEnd);
        y = y(itStart:itEnd);
    end
    
    % attempt to estimate Fs
    % (maybe look for x-axis time labels...)
    dt = x(2) - x(1);
    Fs = 1/dt;
    
    NFFT = length(x);

    if( i == 1 )
        prompt = {'Fs:';'NFFT:'};
        dlgtitle = 'WaveFftStruct Inputs';
        numlines = 1;
        def = {num2str(Fs),num2str(NFFT),num2str(0),};
%         def = {'1','2'}
%         answer = inputdlg(prompt,dlgtitle,numlines,def);
        answer = fft_dialog(def);
    end
%     answer
%     Fs = str2num(answer{1});
%     NFFT = str2num(answer{2});
    Fs = answer{1};
    NFFT = answer{2};
    select_window = answer{3};
    twosided = answer{4};
    fSc = answer{5};
    fU = answer{6};
    fOffset = answer{7};
    
% % %     % selecting a window requires a list of windows to be defined
% % %     files = ls('Z:\sawmillkd\MATLAB\Windows\*.m');
% % %     for iwin = 1:size(files,1)
% % %         [path,windows{iwin},ext] = fileparts(files(iwin,:));
% % %     end
% % %     [Selection,ok] = listdlg('PromptString','Available Windows:',...
% % %         'SelectionMode','single','ListString',windows);
    
%     % default window
%     if(ok == 1)
%         select_window = str2func(windows{Selection});
%     else
%         select_window = @rectwin;
%     end
    
%     WaveFftStruct(y,ones(1,length(x)),Fs,NFFT)
    if(twosided == 0)
        tmpwav = WaveObj(y,select_window(length(x)),Fs,NFFT);
    else
        tmpwav = WaveObj(y,select_window(length(x)),Fs,NFFT,'twosided');
    end
    plot(tmpwav,fSc,fU,fOffset)
    
end


    function choice = fft_dialog(defaults)
        
        selectable_freq_units = {'Hz','kHz','MHz','GHz'};
        
%         hfig = gcf;
        origUnits = get(hfig,'Units');
        set(hfig,'Units','pixels');
        centerpos = get(hfig,'Position');
        set(hfig,'Units',origUnits)
        centerpos = centerpos(1:2) + centerpos(3:4)/2;
        dialog_size = [500 300];
        dialog_pos = [centerpos(1)-dialog_size(1)/2 centerpos(2)-dialog_size(2)/2 dialog_size];
        
        
        d = dialog('Position',dialog_pos,'Name','FFT Settings');
        set(d,'Units','normalized');
        
        dy = 0.065;
        ii = 1;
        txt_fs = uicontrol('Parent',d,'Style','text','units','normalized',...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String','Fs:','HorizontalAlignment','left');

        
        ii = 2;
        inp_fs = uicontrol('Parent',d,'Style','edit','units','normalized', ...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String',defaults{1});
        dropdown_fs_units = uicontrol('Parent',d,'Style','popup', 'units','normalized', ...
            'Position',[0.55 1-ii*dy 0.3 dy], ...
            'String',selectable_freq_units);
        
        ii = 3;
        txt_nfft = uicontrol('Parent',d,'Style','text','units','normalized',...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String','NFFT:','HorizontalAlignment','left');
        ii = 4;
        inp_nfft = uicontrol('Parent',d,'Style','edit','units','normalized', ...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String',defaults{2});
        
        ii = 5;
        txt_window = uicontrol('Parent',d,'Style','text','units','normalized',...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String','Window:','HorizontalAlignment','left');
        
        % selecting a window requires a list of windows to be defined
        ii = 6;
        files = ls('Z:\sawmillkd\MATLAB\Windows\*.m');
        % default window is first in "windows" cell array
        windows{1} = 'rectwin';
        for iwin = 1:size(files,1)
            [path,windows{iwin+1},ext] = fileparts(files(iwin,:));
        end
        dropdown_window = uicontrol('Parent',d,'Style','popup', 'units','normalized', ...
            'Position',[0.1 1-ii*dy 0.8 dy], ...
            'String',windows);

        ii = 7.5;
        txt_plot_freq_units = uicontrol('Parent',d,'Style','text','units','normalized',...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String','Plot freq units:','HorizontalAlignment','left');
        dropdown_plot_freq_units = uicontrol('Parent',d,'Style','popup', 'units','normalized', ...
            'Position',[0.55 1-ii*dy 0.3 dy], ...
            'String',selectable_freq_units);
        
        ii = 9;
        txt_fOffset = uicontrol('Parent',d,'Style','text','units','normalized',...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String','Freq offset:','HorizontalAlignment','left');
        ii = 10;
        inp_fOffset = uicontrol('Parent',d,'Style','edit','units','normalized', ...
            'Position',[0.1 1-ii*dy 0.4 dy], ...
            'String',defaults{1});
        dropdown_fOffset_units = uicontrol('Parent',d,'Style','popup', 'units','normalized', ...
            'Position',[0.55 1-ii*dy 0.3 dy], ...
            'String',selectable_freq_units);
        
        ii = 11.5;
        checkbox_twosided = uicontrol('Parent',d,'Style','radiobutton','units','normalized', ...
            'Position',[0.1 1-ii*dy 0.8 dy], ...
            'String','Two-sided?');
        
        
        
        btn_ok = uicontrol('Parent',d,'Style','pushbutton','units','normalized', ...
            'Position',[0.1 0.05 0.4 dy],...
            'String','Ok','Callback',@get_inputs_and_close);
        
        % init the output, which gets set during the callback
        % "get_inputs_and_close"
        choice{1} = [];
        
        % wait for dialog to close before running to completion
        uiwait(d);
        
        function get_inputs_and_close(varargin)
            
            % get Fs units
            switch get(dropdown_fs_units,'Value')
                case 1
                    % Hz
                    fs_units_sf = 1;
                case 2
                    % kHz
                    fs_units_sf = 1e3;
                case 3
                    % MHz
                    fs_units_sf = 1e6;
                case 4
                    % GHz
                    fs_units_sf = 1e9;
            end
            
            % get plot frequency units
            %   (converts from Hz)
            switch get(dropdown_plot_freq_units,'Value')
                case 1
                    % Hz
                    plot_freq_units_sf = 1;
                case 2
                    % kHz
                    plot_freq_units_sf = 1e-3;
                case 3
                    % MHz
                    plot_freq_units_sf = 1e-6;
                case 4
                    % GHz
                    plot_freq_units_sf = 1e-9;
            end
            plot_freq_units_str = selectable_freq_units{get(dropdown_plot_freq_units,'Value')};
            
            switch get(dropdown_fOffset_units,'Value')
                case 1
                    % Hz
                    fOffset_units_sf = 1;
                case 2
                    % kHz
                    fOffset_units_sf = 1e3;
                case 3
                    % MHz
                    fOffset_units_sf = 1e6;
                case 4
                    % GHz
                    fOffset_units_sf = 1e9;
            end
            
            % get Fs input
            choice{1} = fs_units_sf * str2num(get(inp_fs,'String'));
            choice{2} = str2num(get(inp_nfft,'String'));
            choice{3} = str2func(windows{get(dropdown_window,'Value')});
            choice{4} = get(checkbox_twosided,'Value');
            choice{5} = plot_freq_units_sf;
            choice{6} = plot_freq_units_str;
            choice{7} = plot_freq_units_sf * fOffset_units_sf * str2num(get(inp_fOffset,'String'));
            
            delete(gcf)
        end
        
    end



end
