function catstr = cellstr_cat(st) ...

% concatenate strings from a cell array of strings with spaces inserted
% between the strings
if(iscell(st))
catstr = [];
lengthst = length(st);
for i = 1:lengthst
    catstr = [catstr st{i}];
    if(i < lengthst)
        catstr = [catstr ' '];
    end
end
end

end