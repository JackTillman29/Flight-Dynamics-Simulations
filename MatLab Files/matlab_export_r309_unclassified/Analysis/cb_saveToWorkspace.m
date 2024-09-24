function cb_saveToWorkspace(varargin)

% get handles to objects in axes
hchild = get(gca,'Children');

% 
for k = 1:length(hchild)
    if(length(hchild) > 1)
        x{k} = get(hchild(k),'xdata');
        y{k} = get(hchild(k),'ydata');
    elseif(length(hchild) == 1)
        x = get(hchild(k),'xdata');
        y = get(hchild(k),'ydata');
    end
end

checkLabels = {'Save "x" data to variable named:' ...
               'Save "y" data to variable named:'}; 
varNames = {'x','y'};
items = {x,y};
export2wsdlg(checkLabels,varNames,items,...
             'Save Figure Data to Workspace');


end