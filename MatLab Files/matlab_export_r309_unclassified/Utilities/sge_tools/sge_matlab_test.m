close all; clear all; clc

tool_dir = '/home/holeja/MATLAB/Utilities/sge_tools';
addpath(tool_dir)

%==========================================================================
% create sge manager
sge = sge_matlab;

%==========================================================================
% 1. initialize study
%    ** configures job_name and job_date, 
%       as well as sge directories (out,err)
study_name = 'sge_matlab_test';
sge = sge.initStudy(study_name);


%==========================================================================
% 2. set matlab script to be run on the grid
sge = sge.setMatlabScript('sge_matlab_test_GRIDSCRIPT');


%==========================================================================
% 3. assign matlab results output directory
matlab_results_dir = [tool_dir '/results'];
sge = sge.setMatlabResultsDir(matlab_results_dir);


%==========================================================================
% 4. assign matlab exe shell command options
sge = sge.setMatlabOptions('-nodesktop','-singleCompThread','-sd',tool_dir);


% THE NEXT LINES WOULD BE CALLED MULTIPLE TIMES VIA A LOOP OR SOME SEQUENCE
% OF TRADE PARAMETERS
% 5. assign variables for use in matlab script
numJobs = 20;
for xmean = linspace(0,10,numJobs)
    
    inputs.xmean = xmean;
    inputs.xstd  = 1;
    inputs.N     = 1e4;
    inputs.nbins = 100;

    sge = sge.ReadyJob('inputs',inputs);

    % 6. this job is ready to be submitted to the grid
    sge = sge.SubmitJob;
    
end


% 6. wait for grid runs to finish
while(sge.GridStatus == 0)
%     sge = sge.fetchData;
    pause(0.5)
    clc
    % inspect a particular job
    disp(sge.job(1).matfile)
end

% 7. get data from each grid run
% sge = sge.fetchData;
disp(sge.job(1))



% inputs



