function res = extract_scpi_value(parsed_response,name,offset)
% This function takes the parsed SCPI response and searches for a
% specific field that you specify. It then returns the next value if no
% offset is provided. If you provide an offset, it will use that instead

if(~exist('offset'))
    offset = 1;
end
for k = 1 : length(parsed_response)
    if(strcmpi(name,parsed_response{k}))
        res = parsed_response{k+offset};
        return;
    end
end


% if you got here, something went wrong. Build a useful debug
% message

disp(['SCPI Response']);
for k = 1 : length(parsed_response)
    disp([num2str(k) '      ' parsed_response{k}]);
end
error(['Unable to find key ' name ' in SCPI response']);

end