function display( FLTB )
    disp('FIR Filter Bank Object');
    disp('#           F1           F2');
    for k = 1 : length(FLTB.filters)
        fprintf('%d %12.0f %12.0f\n',k, ...
            get(FLTB.filters{k},'fStartHz'), ...
            get(FLTB.filters{k},'fStopHz') );
    end
end