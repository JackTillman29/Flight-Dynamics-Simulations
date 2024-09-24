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
study_name = 'multitone_mean_to_peak_ratio';
sge = sge.initStudy(study_name);


%==========================================================================
% 2. set matlab script to be run on the grid
sge = sge.setMatlabScript('rand_phase_peak_to_avg_ratio_SCRIPT');


%==========================================================================
% 3. assign matlab results output directory
matlab_results_dir = [pwd '/results'];
sge = sge.setMatlabResultsDir(matlab_results_dir);


%==========================================================================
% 4. assign matlab exe shell command options
sge = sge.setMatlabOptions('-nodesktop','-singleCompThread','-sd',pwd);


% THE NEXT LINES WOULD BE CALLED MULTIPLE TIMES VIA A LOOP OR SOME SEQUENCE
% OF TRADE PARAMETERS
% 5. assign variables for use in matlab script

numJobs = 25;

disp(['# of jobs: ',num2str(numJobs)])

for ijob = 1:numJobs
    
    inputs.nt   = 4;
    inputs.f0   = 1e6;
    inputs.df   = 500;
    inputs.nrep = 1e6;
    
    sge = sge.ReadyJob('inputs',inputs);

    % 6. this job is ready to be submitted to the grid
    sge = sge.SubmitJob;
    
end


% 6. wait for grid runs to finish
tic;
while(sge.GridStatus == 0)
%     sge = sge.fetchData;
    pause(0.5)
    clc
    toc
    % inspect a particular job
%     disp(sge.job(1).matfile)
end


% 7. get data from each grid run
% sge = sge.fetchData;
disp(sge.job(1))

% 8. get matrix of data from all the worker jobs
disp('getting data...')
peakval = sge.CatData('peakval');
meanval = sge.CatData('meanval');
mpr = sge.CatData('mpr');
disp('...done!')



nbins = 10e3;
[hpeak,epeak] = histcounts(peakval,nbins);
[hmean,emean] = histcounts(meanval,nbins);
[hmpr,empr] = histcounts(mpr,nbins);

cpeak = epeak(2:end) - diff(epeak);
cmean = emean(2:end) - diff(emean);
cmpr = empr(2:end) - diff(empr);

figure;
% subplot(2,1,1)
% plot(mpr,'.')
% subplot(2,1,2)
plot(cmpr,hmpr)
xlabel('Mean-to-Peak Ratio')



figure;
plot(10*log10(cmpr),hmpr)
xlabel('Mean-to-Peak Ratio [dB]')

figure;
plot(10*log10(cpeak),hpeak)
hold on;
plot(10*log10(cmean),hmean)
plot(10*log10(cmpr),hmpr)
legend('peak','mean','mpr')



