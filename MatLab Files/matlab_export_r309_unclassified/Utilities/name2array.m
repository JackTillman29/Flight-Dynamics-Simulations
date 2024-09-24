function y = name2array(data,varargin)
    names = {};
    sizes = [];
    for k = 1 : 2 : length(varargin)
        names = [names varargin{k}];
        sizes = [sizes varargin{k+1}];
    end
    
    % Check if sizes agree before proceeding
    if(size(data,2) ~= sum(sizes))
        error(['Size mismatch on input! Data: ' num2str(size(data,2)) ' vs Specified Size: ' num2str(sum(sizes))]);
    end
    
    % Convert sizes to column indices
    cx = vec2idx(sizes);
    
    % populate data
    for k = 1 : length(sizes)
        eval([ 'y.' names{k} ' = data(:,' num2str(cx(k,1)) ':' num2str(cx(k,2)) ');' ]);
    end
end