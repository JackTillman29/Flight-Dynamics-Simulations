function res = check_scpi_value(obj,request,value)
    scpi_sequence(obj,request);
    pause(0.01);
    res = extract_scpi_value( ...
        parse_scpi_response(fread(obj)), value );
    
end