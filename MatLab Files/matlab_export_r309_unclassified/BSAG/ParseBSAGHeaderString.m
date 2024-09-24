function value = ParseBSAGHeaderString(s,param)
    [match_strings,match_tokens]=regexp(s,[param ' = "(\S+)"'],'match','tokens');
    value = str2num(match_tokens{1}{1});
    if(isempty(value))
        value = match_tokens{1}{1};
    end
end