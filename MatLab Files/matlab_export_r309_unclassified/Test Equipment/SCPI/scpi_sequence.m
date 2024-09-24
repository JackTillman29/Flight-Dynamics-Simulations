function scpi_sequence(obj,varargin)
% Syntax: scpi_sequence(obj,'string1','string2',....)
    for k = 1 : length(varargin)
        fprintf(obj,varargin{k});
    end
end