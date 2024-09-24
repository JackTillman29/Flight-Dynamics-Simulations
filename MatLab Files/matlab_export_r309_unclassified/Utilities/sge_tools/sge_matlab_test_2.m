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
study_name = 'sge_matlab_test_2';
sge = sge.initStudy(study_name);


%==========================================================================
% 2. set matlab script to be run on the grid
sge = sge.setMatlabScript('sge_matlab_test_GRIDSCRIPT_2');


%==========================================================================
% 3. assign matlab results output directory
matlab_results_dir = [tool_dir '/results_2'];
sge = sge.setMatlabResultsDir(matlab_results_dir);


%==========================================================================
% 4. assign matlab exe shell command options
sge = sge.setMatlabOptions('-nodesktop','-singleCompThread','-sd',tool_dir);


% THE NEXT LINES WOULD BE CALLED MULTIPLE TIMES VIA A LOOP OR SOME SEQUENCE
% OF TRADE PARAMETERS
% 5. assign variables for use in matlab script

knob_avgPwr_1_dBm_arr = [-5:1:0];
knob_avgPwr_2_dBm_arr = [-5:1:0];

n1 = length(knob_avgPwr_1_dBm_arr);
n2 = length(knob_avgPwr_2_dBm_arr);

numJobs = n1 * n2;

disp(['# of jobs: ',num2str(numJobs)])

k1 = 0;
for knob_avgPwr_1_dBm = knob_avgPwr_1_dBm_arr
    k1 = k1 + 1;
    k2 = 0;
for knob_avgPwr_2_dBm = knob_avgPwr_1_dBm_arr
    k2 = k2 + 1;
    
    inputs.Rout  = 50;
    inputs.N     = 1e6;
    inputs.nreps = 5;
    
    inputs.knob_avgPwr_1_dBm = knob_avgPwr_1_dBm;
    inputs.knob_avgPwr_2_dBm = knob_avgPwr_2_dBm;
    
    inputs.k1 = k1;
    inputs.k2 = k2;
    inputs.n1 = n1;
    inputs.n2 = n2;
    inputs.meanPwr_dBm = 'eval.zeros(n1,n2)';
    
    sge = sge.ReadyJob('inputs',inputs);

    % 6. this job is ready to be submitted to the grid
    sge = sge.SubmitJob;
    
end
end


% 6. wait for grid runs to finish
tic;
while(sge.GridStatus == 0)
%     sge = sge.fetchData;
    pause(0.5)
    clc
    toc
    % inspect a particular job
    disp(sge.job(1).matfile)
end

% 7. get data from each grid run
% sge = sge.fetchData;
disp(sge.job(1))

% 8. get matrix of data from all the worker jobs
meanPwr_dBm = sge.GetData('meanPwr_dBm');

figure; imagesc(knob_avgPwr_1_dBm_arr,knob_avgPwr_2_dBm_arr,meanPwr_dBm);
set(gca,'YDir','normal')
colorbar
xlabel('Power of Signal 2 [dBm]')
ylabel('Power of Signal 1 [dBm]')
title('Mean Power of two summed complex AWGN signals')

% inputs



