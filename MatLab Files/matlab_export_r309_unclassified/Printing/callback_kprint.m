function callback_kprint(varargin)

if(nargin > 0)
%     for k = 1:2:length(varargin)
%         thisString = varargin{k};
%         thisValue  = varargin{k+1};
%         
%         switch thisString
%             case 'I Feel Lucky'
                
%     end

    optionString = varargin{1};

else
    optionString = 'Dialog';

end
    
hfigchild = get(gcf,'Children');

% default to figure X for filename
try
    filename = ['Figure ' num2str(get(gcf,'Number'))];
catch
    filename = ['Figure ' num2str(gcf)];
end

% find figure's axes object(s) and search for a title string
try
    for i = length(hfigchild):-1:1
        child = hfigchild(i);
        if(~isempty(strfind(child.Type,'axes')))
            if(~isempty(child.Title.String))
                % if multiple strings are in the title, then prompt user to
                % select which string to use as filename
                if(iscell(child.Title.String))
                    titlestr = child.Title.String;
%                     titlestr{end+1} = [titlestr{:}];
                    allTitles = cellstr_cat(titlestr);
                    titlestr{end+1} = allTitles;
                    titlestr{end+1} = ['Figure ' num2str(get(gcf,'Number'))];
                    
                    % "optionString" is an optional input argument to this function
                    % that defaults to 'Dialog'
                    switch optionString
                        case 'Dialog'
                            idx = listdlg('PromptString','Select a filename:',...
                                'SelectionMode','single','ListString',titlestr,...
                                'ListSize',[400 150]);
                            filename = [titlestr{idx}];
                        case 'I Feel Lucky'
                            filename = [allTitles];
                        otherwise
                            disp('Unrecognized optional input argument, "optionString" in ''callback_kprint''')
                    end
                else
                    filename = [child.Title.String];
                end
            end
        end      
    end
catch
    for i = length(hfigchild):-1:1
    child = hfigchild(i);
    if(~isempty(strfind(get(child,'Type'),'axes')))
        titlestr = get(get(child,'Title'),'String');
        if(~isempty(titlestr))
            % if multiple strings are in the title, then prompt user to
            % select which string to use as filename
            if(iscell(titlestr))
%                 titlestr{end+1} = [titlestr{:}];
                allTitles = cellstr_cat(titlestr);
                titlestr{end+1} = allTitles;
                titlestr{end+1} = ['Figure ' num2str(gcf)];
                switch optionString
                    case 'Dialog'
                        idx = listdlg('PromptString','Select a filename:',...
                                'SelectionMode','single','ListString',titlestr,...
                                'ListSize',[400 150]);
                        if(exist('idx','var'))
                            filename = [titlestr{idx}];
                        else
                            filename = [titlestr];
                        end
                    case 'I Feel Lucky'
                        filename = [allTitles];
                    otherwise
                        disp('Unrecognized optional input argument, "optionString" in ''callback_kprint''')
                end
            end
            
            break
        end
    end      
    end
    
    
end

try
    filename = [filename datestr(datetime)];
catch
    filename = [filename datestr(clock)];
end

filename=strrep(filename,'/','_');
filename=strrep(filename,'\','_');
filename=strrep(filename,':',' ');


% define the default fileExt here.
fileExt = '.png';

switch optionString
    case 'Dialog'
        a=inputdlg({'Image Name','DPI'},'Options',1,{[filename fileExt],'180'});
    case 'I Feel Lucky'
        % no dialog boxes, simply efficient.
        a{1} = [filename fileExt];
        a{2} = '180'; % DPI
    otherwise
        disp('Unrecognized optional input argument, "optionString" in ''callback_kprint''')
end
        

try
    kprint(get(gcf,'Number'),a{1},str2num(a{2}));
catch
    kprint(gcf,a{1},str2num(a{2}));
end

end