function cb_modify_colorbar_caxis(varargin)

hfig = varargin{1};
hKeyData = varargin{2};

hcb = findobj('Parent',hfig,'Type','colorbar');
hax = findobj('Parent',hfig,'Type','axes');
oldflag = 0;
if(isempty(hcb))
    % for backwards compat.
    hcb = findobj('Parent',hfig,'Type','axes','Tag','Colorbar');
    hax = findobj('Parent',hfig,'Type','axes','-not','Tag','Colorbar');
    oldflag = 1;
end



% ensure that the colorbar exists
if(~isempty(hcb))
    if(oldflag == 0)
        clims = get(hcb,'Limits');
    else
        clims = get(hcb,'Ylim');
    end

try
    if(strcmp(get(hKeyData,'Modifier'),{'shift' 'control'}) & [1 1])
        % increase upper caxis limit
        if(strcmp(get(hKeyData,'Key'),'uparrow'))
            newclims = [clims(1) ceil(clims(2)+1)];
        % decrease upper caxis limit
        elseif(strcmp(get(hKeyData,'Key'),'downarrow'))
            newclims = [clims(1) ceil(clims(2)-1)];
        end
           
    elseif(strcmp(get(hKeyData,'Modifier'),'control'))
        % increase lower caxis limit
        if(strcmp(get(hKeyData,'Key'),'uparrow'))
            newclims = [floor(clims(1)+1) clims(2)];
        % decrease lower caxis limit
        elseif(strcmp(get(hKeyData,'Key'),'downarrow'))
            newclims = [floor(clims(1)-1) clims(2)];
        end 
    end
catch
    error('issue with parsing KeyData object...')
end

% prevent the upper limits from being less than or equal to the lower limit
if(exist('newclims'))
    if(newclims(2) <= newclims(1))
        caxis(clims)
    else
        caxis(newclims)
    end
end

end


% end
