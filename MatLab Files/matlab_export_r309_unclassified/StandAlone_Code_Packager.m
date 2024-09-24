function sa_package(file,output_folder)
    % Syntax: StandAlone_Code_Packager('input_file.m','output_folder')
    warning(['Make sure that your code runs before calling this function '...
        'to ensure all needed functions are on the MATLAB search path. ' ...
        'Otherwise, files may be missing.'])
    
    dep = matlab.codetools.requiredFilesAndProducts(file);
    mkdir(output_folder);
    [fin_path,fin_name,fin_ext]=fileparts(file);
    for k = 1 : length(dep)
        [d_path,d_name,d_ext]=fileparts(dep{k});
        
        % Check if there is a local copy of the dependency file in the
        % input folder
        if(exist([fin_path '\' d_name '.m'],'file'))
            disp(['Copying ' fin_path '\' d_name '.m']);
            copyfile([fin_path '\' d_name '.m'],output_folder);
        else
            disp(['Copying ' dep{k}]);
            copyfile(dep{k},output_folder);
        end
    end
end