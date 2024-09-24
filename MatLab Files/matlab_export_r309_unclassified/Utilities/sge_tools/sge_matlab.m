classdef sge_matlab
    
    properties
        
        % user inputs that might be nice to change
        save_matfile = 1; % 1: yes, 0: no
        
        % define executable paths
        matlab_exe = '';
        shell_exe  = '';
        sge_exe    = '';
        
        study_name = '';
        study_date = '';
        
        % define matlab script to run on grid
        matlab_script = '';
        
        user           = '';
        home_dir       = '';
        sge_out_dir    = '';
        sge_err_dir    = '';
        matlab_results_dir = ''; % where any data is saved
        
        matlab_options     = {}; % cell array of strings
        % THIS VARIABLE IS HOW INPUTS ARE PASSED INTO THE WORKER MATLAB
        matlab_run_options = {}; % cell array of strings
        
        
        matlab_options_str        = ''; % semi-private (let class set)
        matlab_run_options_str    = ''; % semi-private (let class set)
        matlab_full_shell_command = ''; % semi-private (let class set)
        sge_job_script            = ''; % semi-private (let class set)
        sge_full_shell_command    = ''; % semi-private (let class set)
        sge_other_opts            = '';
        
        job_output = [];
        job_count = 0;
        job = struct([]); %struct('number','matfilename','matfile');
        
        run_option_eval_prefix_regexp = 'eval\.'; % REGEXP EXPRESSION
    end
    
    
    %% MAIN INITIALIZATION FUNCTIONS
    %======================================================================
    % MAIN INITIALIZATION FUNCTIONS
    %======================================================================
    methods
        
        
        function this = sge_matlab(varargin)
            % ### CONSTRUCTOR ###
            %   varargin{1}: a *.m file that contains all the defaults for
            %                use by this class object, defining:
            %                - matlab_exe
            %                - sge_exe
            %                - shell_exe
            % how to use this class:
            %    1. Initialize study:
            %           sge = sge.initStudy('give_me_a_name');
            %    2. set matlab script to be run on the grid
            %           sge = sge.setMatlabScript('any_std_matlab_script');
            %    3. assign matlab results output directory
            %           sge = sge.setMatlabResultsDir('path/to/dir/');
            %    4. assign matlab exe shell command options
            %           sge = sge.setMatlabOptions(opt1,opt2);
            %           sge = sge.setMatlabOptions({opt1,opt2});
            if(length(varargin) > 0)
                defaults_file = varargin{1};
            end
            if(exist('defaults_file'))
                % defaults_file should be a MATLAB script that sets any/all
                % of the variables IN THIS SECTION.
                run(defaults_file)
            else
                % defaults for location of this local copy
                this.matlab_exe = '/usr/local/bin/matlab';
                this.sge_exe    = 'qsub';
                this.shell_exe  = '#/bin/bash';
            end
            
            [junk,username] = unix('echo $USER');
            [junk,home_dir] = unix('echo $HOME');
            
            this.home_dir = home_dir(1:end-1);
            this.user     = username(1:end-1);
            
        end
        
        function this = setMatlabPath(this,matlab_exe)
            if(isa(matlab_exe,'char') || isa(matlab_exe,'string'))
                this.matlab_exe = matlab_exe;
            else
                error(['matlab_exe should be a character array'])
            end
        end
        
        function this = setMatlabScript(this,matlab_script)
            if(isa(matlab_script,'char') || isa(matlab_script,'string'))
                % check for extension (at end of char array)
                if(isempty(regexpi(matlab_script,'.m$')))
                    matlab_script = [matlab_script '.m'];
                end
                % check to see if file exists
                %    2: is a file with extension .m, .mlx, or .mlapp, or
                %       name is the name of a file with a non-registered
                %       file extension (.mat, .fig, .txt).
                if(exist(matlab_script) == 2)
                    this.matlab_script = matlab_script;
                else
                    error([matlab_script,...
                        ' does not exist. check filename and/or path.'])
                end
            else
                error(['matlab_script should be a character array'])
            end
        end
        
        function this = initStudy(this,study_name)
            % this = initStudy(this,study_name)
            %    study_name: char array;
            %                must not contain any special characters, and
            %                whitespace will be replaced with underlines.
            %                "study_name" will be used by this class to
            %                build directories and scripts, so this name
            %                must not use any special characters reserved
            %                by the OS. Any whitespace will be replaced
            %                with underlines so that the filenames of the
            %                scripts built with "study_name" can be
            %                executed by the OS
            % assigns study_name and uses this along with a date/time stamp
            % to set the grid output and error directories, as well as the
            % sge jobe script name.
            if(isa(study_name,'char') || isa(study_name,'string'))
                this.study_name = strrep(study_name,' ','_');
                this.study_date = datestr(now,'yyyymmdd_HHMMSS.FFF');
                this = this.setGridOutputDir(...
                    [this.home_dir '/grido_',this.study_name]);
                this = this.setGridErrorDir(...
                    [this.home_dir '/gride_',this.study_name]);
                this.sge_job_script = [this.study_name '.bash'];
            else
                error(['matlab_script should be a character array'])
            end
        end
        
        function this = setMatlabResultsDir(this,matlab_results_dir,varargin)
            if(this.save_matfile)
                if(isa(matlab_results_dir,'char') ...
                        || isa(matlab_results_dir,'string'))
                    this.matlab_results_dir = matlab_results_dir;
                    % make directory (if doesn't exist)
                    mkdir(this.matlab_results_dir)
                    % clear directory of all existing content
                    %delete([this.matlab_results_dir,'/*'])
                else
                    error(['matlab_results_dir should be a character array'])
                end
            else
                error('User requested NOT to save matfile ("save_matfile" field is set to zero)')
            end
        end
        
        function this = setGridOutputDir(this,sge_out_dir)
            if(isa(sge_out_dir,'char') || isa(sge_out_dir,'string'))
                this.sge_out_dir = sge_out_dir;
                % make directory (if doesn't exist)
                mkdir(this.sge_out_dir)
                % clear directory of all existing content
                delete([this.sge_out_dir,'/*'])
            else
                error(['sge_out_dir should be a character array'])
            end
        end
        
        function this = setGridErrorDir(this,sge_err_dir)
            if(isa(sge_err_dir,'char') || isa(sge_err_dir,'string'))
                this.sge_err_dir = sge_err_dir;
                % make directory (if doesn't exist)
                mkdir(this.sge_err_dir)
                % clear directory of all existing content
                delete([this.sge_err_dir,'/*'])
            else
                error(['sge_err_dir should be a character array'])
            end
        end
        
        function this = setShellCommand(this,shell_exe)
            if(isa(shell_exe,'char') || isa(shell_exe,'string'))
                this.shell_exe = shell_exe;
            else
                error(['shell_exe should be a character array'])
            end
        end
        
        function this = setSgeShellCommand(this,sge_exe)
            if(isa(sge_exe,'char') || isa(sge_exe,'string'))
                this.sge_exe = sge_exe;
            else
                error(['sge_exe should be a character array'])
            end
        end
        
        function this = setMatlabOptions(this,varargin)
            % this = SetOptionsString(this,options)
            % options needs to be a CELL array of STRINGS, each string
            % containing the option flag and/or user inputs
            % EXAMPLES:
            %   matlab_options= {'-nodesktop'}
            %   matlab_options= {'-nodesktop','-sd',
            %                   '"/path/to/starting/directory/"'}
            %     ** for options expecting user input, you can put the user
            %     input in the same string as the option flag, or you can
            %     put it in the next cell element (either should work):
            %
            %     matlab_options= {'-sd','"/path/to/starting/directory/"'}
            %         OR
            %     matlab_options= {'-sd "/path/to/starting/directory/"'}
            if(isa(varargin,'cell'))
                this.matlab_options     = varargin;
                this.matlab_options_str = ...
                    BuildOptionsString(this,varargin);
            else
                error(['matlab_options need to be input as a ',...
                    'cell array of strings. ( ',...
                    'examples: {''-nodesktop'', ''-sd'',',...
                    ' ''"/path/to/starting/dir/"''} )'])
            end
        end
        
    end
    
    %% DATA PASSING FUNCTIONS (MASTER --> WORKER)
    %======================================================================
    % DATA PASSING FUNCTIONS (FROM MASTER MATLAB TO WORKER MATLABS)
    %======================================================================
    methods
        
        % MAIN DISCUSSION =================================================
        % Since the purpose of this class is to run jobs independently in
        % parallel on a multicore machine on multiple instantiations of
        % MATLAB, the need for a "clear all" line at the top of our scripts
        % is not necessary because each MATLAB instance starts FRESH and
        % only accepts values you tell it from the RUN OPTION command "-r
        % ...". For ease of passing variables INTO this class, it is
        % recommended that any variables that you want in the WORKER MATLAB
        % be instantiated and defined in a set of STRUCTs that are passed
        % into the class. The class will unpack the STRUCTs and create
        % MATLAB commands for the RUN OPTION that will be executed on the
        % WORKER MATLAB to create these variables and their values.
        %
        %                   (pass into)    (strings)
        %  STRUCT myInputs --> class --|-->var1 = ...; --|
        %                              |-->var2 = ...; --|----->
        %                              |-->var3 = ...; --|
        %    --> concat string --> run_option_combined_string
        %
        %    matlab -r "run_option_combined_string"
        %                |
        %                V
        %           where "run_option_combined_string" is actually expanded
        %           out in the shell script command.
        
        
        % Any values that you want in the worker matlab from the master
        % matlab need to be passed into the '-r "......"' option
        % explicitly, one by one. One solution is to build a "gridinputs"
        % structure containing all the variables you want in the worker and
        % then have this function build each "<prev cmd>; thisvar =
        % values;" and append to the -r option string for the shell script
        % ** If you want to pass a large array, want to build matrices, or
        % compute things which depend on previously defined variables in
        % the structure IN THE RUN STRING (which is RUN by WORKER MATLAB),
        % then you can actually define a CHARACTER ARRAY with the PREFIX
        % "eval\." (a REGEXP expression) to tell the function below that
        % you want the string to go INSIDE an eval() command.
        
        % SIDEBAR =========================================================
        % If the script you are running ON THE GRID does not need to have a
        % "clear all;" or some variant at the top, then you can simply let
        % the variables
        
        function this = BuildRunOptionString(this,varargin)
            % varargin are STRING-STRUCT pairs where the STRING (in
            % varargin{k}) is the name STRUCT variable (in varargin{k+1})
            this.matlab_run_options = cell(0); % init to size 0
            for ivar = 1:2:length(varargin)
                if(isa(varargin{ivar},'char') && ...
                        isa(varargin{ivar+1},'struct'))
                    %   need to set each field of varargin{ivar+1} in run
                    %   option string
                    fields = fieldnames(varargin{ivar+1});
                    for k = 1:length(fields)
                        fielddat = varargin{ivar+1}.(fields{k});
                        if(isscalar(fielddat))
                            % field content is SCALAR
                            this.matlab_run_options{k} = ...
                                [fields{k},' = ',num2str(fielddat),';'];
                        elseif(isa(fielddat,'char'))
                            % field content is CHAR string
                            [i1,i2] = regexpi(fielddat,...
                                ['^',this.run_option_eval_prefix_regexp]);
                            if(~isempty(i1))
                                % create eval statement
                                %   (regexp run_option_eval_prefix found)
                                this.matlab_run_options{k} = ...
                                    ['eval([''',fields{k},' = ', ...
                                        fielddat((i2+1):end),''']);'];
                            else
                                % set variable to be a CHAR string
                                this.matlab_run_options{k} = ...
                                    [fields{k},' = ''',fielddat,''';'];
                            end
                        elseif(isa(fielddat,'string'))
                            % field content is STRING type
                            %  (string type introduced in R2016b, may have
                            %  compatibility issues here...)
                            this.matlab_run_options{k} = ...
                                [fields{k},' = "',fielddat,'";'];
                        elseif(isnumeric(fielddat) & (size(fielddat,1) == 1) & (size(fielddat,2) > 1))
                            % field content is a ROW VECTOR
                            this.matlab_run_options{k} = ...
                                [fields{k},' = [',num2str(fielddat),'];'];
                        elseif(isnumeric(fielddat) & (size(fielddat,1) > 1) & (size(fielddat,2) == 1))
                            % field content is a COL VECTOR
                            this.matlab_run_options{k} = ...
                                [fields{k},' = [',num2str(fielddat.'),'].'];
                        end
                    end
                    
                else
                    error(['all inputs need to be STRUCTs'])
                end
            end
            % now that it's finished, build string from run options
            this.matlab_run_options_str = ...
                this.BuildOptionsString(this.matlab_run_options);
        end
        
        
        function this = BuildMatlabFullShellCommand(this)
            %
            % This command assumes the following class variables have been
            % defined:
            %  matlab_exe
            %  matlab_script
            %  matlab_options
            %  matlab_run_options
            
            run_mfile_cmd = ['run(''',this.matlab_script,'''); '];
            
            if(this.save_matfile == 1)
                save_matfile_cmd = [...
                    'save([''',this.job(this.job_count).matfilename,''']); '];
            else
                save_matfile_cmd = '';
            end
            %save_matfile_cmd = '';
            
            % init random number generator
            init_rng_cmd = ['rng(',num2str(mod(datenum(now),1)*1e8,'%50.30f'),'); '];
            
            % build matlab command for bash script
            this.matlab_full_shell_command = ...
                [this.matlab_exe this.matlab_options_str, ...
                ' -r "', ...
                this.matlab_run_options_str,...
                init_rng_cmd, ...
                run_mfile_cmd, ...
                save_matfile_cmd, ...
                'quit;"'];
            % TODO: this.matlab_script could be absorbed into the
            % "this.run_option" definitions as the input structures are
            % built up using the EVAL paradigm described above at the top
            % of this section.
        end
        
    end
    
    %% QUERY GRID FUNCTIONS
    %======================================================================
    % QUERY GRID FUNCTIONS
    %======================================================================
    methods
        function status = GridStatus(this,option)
            if(~exist('option'))
                option = 'summary';
            end
            
            ipend = 0;
            irun  = 0;
            idone = 0;
            
            for ijob = 1:length(this.job)
                status = JobStatus(this, ijob, option);
                switch status
                    case 'pending'
                        ipend = ipend + 1;
                    case 'running'
                        irun = irun + 1;
                    case 'done'
                        idone = idone + 1;
                end
            end
            
            disp(['Number of jobs: ',num2str(this.job_count)])
            disp(['Jobs running:   ',num2str(irun)])
            disp(['Jobs waiting:   ',num2str(ipend)])
            disp(['Jobs finished:  ',num2str(idone)])
            
            if(idone < this.job_count)
                status = 0;
            else
                status = 1;
            end
            
%             switch lower(option)
%                 case {'summary','s'}
%                     % METHOD 1 - LOOP THROUGH JOBS ONE BY ONE
% %                     job_done = zeros(1,this.job_count);
% %                     for k = 1:this.job_count
% %                         [job_done(k),cmdout]=unix([...
% %                             'qstat -j ',num2str(jobnum(k))]);
% %                     end
%                     % METHOD 2 - USE CLI TOOLS TO GET LIST OF JOBS ALL AT
%                     %            ONCE (filtering out junk lines)
%                     [status,jobs_run] = unix([...
%                         'qstat -s r -u $USER | ',...
%                         'awk ''/[0-9]/{print $1}'' | wc -l']);
%                     [status,jobs_wait] = unix([...
%                         'qstat -s p -u $USER | ',...
%                         'awk ''/[0-9]/{print $1}'' | wc -l']);
%                     jobs_run  = str2num(jobs_run);
%                     jobs_wait = str2num(jobs_wait);
%                     jobs_finished = this.job_count - jobs_run - jobs_wait;
%                     disp(['Number of jobs: ',num2str(this.job_count)])
%                     disp(['Jobs running:   ',num2str(jobs_run)])
%                     disp(['Jobs waiting:   ',num2str(jobs_wait)])
%                     disp(['Jobs finished:  ',num2str(jobs_finished)])
%                     
%                     if(jobs_finished < this.job_count)
%                         status = 0;
%                     else
%                         status = 1;
%                     end
%                     
%                 case {'detailed','details','d'}
%                     for k = 1:this.job_count
%                         [status,cmdout] = unix([...
%                             'qstat -j ',num2str(this.job(k).sge_job_id),...
%                             ' | awk ''/usage/{x=$3; gsub("cpu=","",x);',...
%                             'gsub(",","",x);}''']);
%                         if(status == 1)
%                             if(isempty(cmdout))
%                                 cpu_time = 'pending';
%                             else
%                                 cpu_time = cmdout;
%                             end
%                             disp(['job id:  ',num2str(k),'  ',cpu_time,])
%                             disp(['sge job id:  ',num2str(this.job(k).sge_job_id),'  ',cpu_time])
%                         else
%                             disp(['status returned 0 for sge job # ',num2str(this.job(k).sge_job_id)])
%                         end
%                     end
%                 otherwise
%                     warning(['unknown option passed in. ',...
%                         'exiting without doing anything'])
%             end
                    
        end
        
        function status = JobStatus(this, ijob, option)
            if(~exist('option'))
                option = 'summary';
            end
            switch lower(option)
                case {'summary','s'}
                    % METHOD 1 - LOOP THROUGH JOBS ONE BY ONE
%                     job_done = zeros(1,this.job_count);
%                     for k = 1:this.job_count
%                         [job_done(k),cmdout]=unix([...
%                             'qstat -j ',num2str(jobnum(k))]);
%                     end
                    % METHOD 2 - USE CLI TOOLS TO GET LIST OF JOBS ALL AT
                    %            ONCE (filtering out junk lines)
                    [status,cmdout] = unix([...
                        'qstat -j ',num2str(this.job(ijob).sge_job_id)]);
                    
                    % look at cmdout for these kinds of strings:
                    %   **NOTE: this will not work in the future if SGE
                    %   changes its output formatting or text.
                    if( strfind(cmdout,'jobs do not exist') ~= 0 )
                        status = 'done';
                    elseif( strfind(cmdout,'cpu=') ~= 0 )
                        status = 'running';
                    else
                        status = 'pending';
                    end
                        
                    
%                     [status,jobs_run] = unix([...
%                         'qstat -s r -u $USER | ',...
%                         'awk ''/[0-9]/{print $1}'' | wc -l']);
%                     [status,jobs_wait] = unix([...
%                         'qstat -s p -u $USER | ',...
%                         'awk ''/[0-9]/{print $1}'' | wc -l']);
%                     jobs_run  = str2num(jobs_run);
%                     jobs_wait = str2num(jobs_wait);
%                     jobs_finished = this.job_count - jobs_run - jobs_wait;
%                     disp(['Number of jobs: ',num2str(this.job_count)])
%                     disp(['Jobs running:   ',num2str(jobs_run)])
%                     disp(['Jobs waiting:   ',num2str(jobs_wait)])
%                     disp(['Jobs finished:  ',num2str(jobs_finished)])
%                     
%                     if(jobs_finished < this.job_count)
%                         status = 0;
%                     else
%                         status = 1;
%                     end
                    
%                 case {'detailed','details','d'}
%                     for k = 1:this.job_count
%                         [status,cmdout] = unix([...
%                             'qstat -j ',num2str(this.job(k).sge_job_id),...
%                             ' | awk ''/usage/{x=$3; gsub("cpu=","",x);',...
%                             'gsub(",","",x);}''']);
%                         if(status == 1)
%                             if(isempty(cmdout))
%                                 cpu_time = 'pending';
%                             else
%                                 cpu_time = cmdout;
%                             end
%                             disp(['job id:  ',num2str(k),'  ',cpu_time,])
%                             disp(['sge job id:  ',num2str(this.job(k).sge_job_id),'  ',cpu_time])
%                         else
%                             disp(['status returned 0 for sge job # ',num2str(this.job(k).sge_job_id)])
%                         end
%                     end
%                 otherwise
%                     warning(['unknown option passed in. ',...
%                         'exiting without doing anything'])
            end
        end
    end
    
    %% DATA PASSING FUNCTIONS (WORKERS --> MASTER)
    %======================================================================
    % DATA COMBINING FUNCTIONS (GET DATA BACK FROM THE WORKERS)
    %======================================================================
    methods
        
        function data = GetData(this,jobnum,varname)
            if(~exist('varname'))
                error('need to pass in a variable name')
            end
            
            try
                tmp = this.job(jobnum).matfile.(varname);
                if(~exist('data'))
                    data = tmp;
                else
                    data = data + tmp;
                end
            catch
                warning(['data from job(',num2str(jobnum),...
                    ') is either still running, encountered ',...
                    'an error, or the desired data doesn''t ',...
                    'exist'])
            end
        end
        
        function data = GetDataAll(this,varname)
            
            if(~exist('varname'))
                error('need to pass in a variable name')
            end
            
            for ijob = 1:this.job_count
                try
                    tmp = this.job(ijob).matfile.(varname);
                    if(~exist('data'))
                        data = tmp;
                    else
                        data = data + tmp;
                    end
                catch
                    warning(['data from job(',num2str(ijob),...
                        ') is either still running, encountered ',...
                        'an error, or the desired data doesn''t ',...
                        'exist'])
                end
            end
        end
        
        function data = CatData(this,varname)
            
            if(~exist('varname'))
                error('need to pass in a variable name')
            end
            
            for ijob = 1:this.job_count
                try
                    tmp = this.job(ijob).matfile.(varname);
                    if(~exist('data'))
                        data = tmp;
                    else
                        % ROW VECTOR
                        if( (size(tmp,2) > 1) & (size(tmp,1) == 1) )
                            data = [data tmp];
                        % COLUMN VECTOR
                        elseif( (size(tmp,2) == 1) & (size(tmp,1) > 1) )
                            data = [data; tmp];
                        % MATRIX DATA (TRY CONCATENATING)
                        elseif( (size(tmp,2) > 1) & (size(tmp,1) > 1) )
                            data = [data tmp];
                        % SCALAR DATA
                        else
                            data = [data tmp];
                        end
                    end
                catch
                    warning(['data from job(',num2str(ijob),...
                        ') is either still running, encountered ',...
                        'an error, or the desired data doesn''t ',...
                        'exist'])
                end
            end
        end
        
    end
    
    %% SHELL SCRIPT CREATION FUNCTIONS
    %======================================================================
    % SHELL SCRIPT CREATION FUNCTIONS
    %======================================================================
    methods
        function CreateJobShellScript(this)
            % builds a temporary script
            fid = fopen(this.sge_job_script,'w');
            fprintf(fid,'%s\n',this.shell_exe);
            fprintf(fid,'%s\n',this.matlab_full_shell_command);
            fclose(fid);
        end
        
        function this = BuildSgeFullShellCommand(this)
            sge_job_cmd = [this.sge_exe,' '];
            if(~isempty(this.sge_err_dir))
                sge_job_cmd = [sge_job_cmd '-e ', this.sge_err_dir,' '];
            end
            if(~isempty(this.sge_out_dir))
                sge_job_cmd = [sge_job_cmd '-o ', this.sge_out_dir,' '];
            end
            if(~isempty(this.sge_other_opts))
                sge_job_cmd = [sge_job_cmd this.sge_other_opts,' '];
            end
            sge_job_cmd = [sge_job_cmd,' ',this.sge_job_script];
            
            this.sge_full_shell_command = sge_job_cmd;
        end
        
        function this = ReadyJob(this,varargin)
            % varargin are STRING-STRUCT pairs where the STRING (in
            % varargin{k}) is the name STRUCT variable (in varargin{k+1})
            
            % increment job counter
            this.job_count = this.job_count + 1;
            
            % initialize job structure
            this.job(this.job_count).number      = this.job_count;
            this.job(this.job_count).matfilename = [...
                    this.matlab_results_dir,'/', ...
                    this.study_name,'_',num2str(this.job_count),'.mat'];
            
            this.job(this.job_count).matfile = matfile( ...
                    this.job(this.job_count).matfilename, ...
                    'Writable',logical(0));
            
            
            % PASS INPUT PARAMETERS (STRING-STRUCT PAIRS) TO BE USED BY JOB
            this = this.BuildRunOptionString(varargin{:});
            
            % USING CURRENT PARAMETERS IN CLASS OBJECT...
            
            % build full matlab shell command
            this = this.BuildMatlabFullShellCommand;
            
            % build job shell script
            this.CreateJobShellScript;
            
            % build sge shell command
            this = this.BuildSgeFullShellCommand;
            
            this.job(this.job_count).matlabFullShellCommand = this.matlab_full_shell_command;
            this.job(this.job_count).fullShellCommand = this.sge_full_shell_command;
            
            pause(0.1)
            
        end
        
        function this = ReadyJobWithId(this,id,varargin)
            
            this.job_count = id - 1;
            
            this = this.ReadyJob(varargin{:});
            
        end
        
        function this = SubmitJob(this)
            
            % submit job to grid through SGE
            [status,cmdout] = unix(this.sge_full_shell_command);
            
            try
            % get job number for this grid submission
            [status,jobnumstr] = unix(...
                ['echo "',cmdout(1:end-1),'" | awk ''//{print $3}''']);
            
            this.job(this.job_count).sge_job_id = str2num(jobnumstr);

            disp(['submitted job (SGE) ',...
                num2str(this.job(this.job_count).sge_job_id),...
                ' (job ',num2str(this.job_count),' - ',this.study_name,')'])
            catch
                disp(this.job_count)
            end
            
        end
        
    end
    
    %% UTILITIES
    %======================================================================
    % UTILITIES
    %======================================================================
    methods
        function optstr = BuildOptionsString(this,matlab_options)
            % optstr = BuildOptionsString(this,options)
            % normally used internally
            % use previously defined "options" if no options are supplied
            if(~exist('matlab_options') & ~isempty(this.matlab_options))
                matlab_options = this.matlab_options;
%             else
%                 warning(['no options supplied... ',...
%                     'recommended to run with "-nodesktop" and make ',...
%                     'explicit a starting directory using ',...
%                     '"-sd /path/to/starting/dir"'])
            end
            
            % INSERT SPACES BETWEEN COMMAND FOR FINAL STRING (NOT CELL
            % ARRAY)
            options2 = cell(1,2*length(matlab_options));
            options2(1:2:end) = {' '}; % add spaces between options
            options2(2:2:end) = matlab_options;
            optstr = cell2mat(options2); % char array output with spaces
                                         % between commands
        end
    end
    
end