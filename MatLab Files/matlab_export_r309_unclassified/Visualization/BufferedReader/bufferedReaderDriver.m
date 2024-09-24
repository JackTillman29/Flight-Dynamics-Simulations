close all;
clear all;
clear classes;

x = bufferedReader( ...
    '<full_path_and_binary_filename>', ...    % filename
    4, ...                      % record header size (bytes)
    4, ...                      % record footer size (bytes)
    3*1000, ...                    % page memory size (# of records)
    'single');                  % precision of data (single or double)



plot(x,@fmt_waveform2cabs);
