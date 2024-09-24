function addpathIncludeSubdirs(topdir,varargin)
% inputs:
%   topdir      -  directory to add to path
%   excludeDirs [optional] - cell array of strings of dirs to not include


if(length(varargin)==1)
    excludeDirs = varargin{1};
elseif(length(varargin)>1)
    error('check inputs')
else
    excludeDirs = [];
end

folders = ls(topdir);
irem = 0;
for k = 1:size(folders,1)
    folder = folders(k,:);
    
    % remove any hidden files/folders from addpath
    i = strfind(folder,'.');
    if(~isempty(i))
        if(i(1)==1)
            irem = irem + 1;
            removeRow(irem) = k;
%             folders(k,:) = [];
        end
    end
    
    % remove scripts from list of folders to add
    i = strfind(folder,'.m');
    if(~isempty(i))
        irem = irem + 1;
        removeRow(irem) = k;
%         folders(k,:) = [];
    end
    
    i = strfind(folder,'old');
    if(~isempty(i))
        if(i(1) == 1)
            irem = irem + 1;
            removeRow(irem) = k;
        end
    end
    
    % remove exclude dirs
    for k = 1:length(excludeDirs)
        i = strfind(folder,excludeDirs(k));
        if(~isempty(i))
            irem = irem + 1;
            removeRow(irem) = k;
        end
    end
        
    
end
folders(removeRow,:) = [];

disp(['top directory: ',topdir])
disp('sub-directories added:')
for k = 1:size(folders,1)
    disp(['   ' folders(k,:)])
    addpath([ topdir strtrim(folders(k,:)) ])
end


end