close all;
fclose all;
clear;
clc;

% Load file
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2193 D20 NBG Baseline V4_1 STIM.dat';
%BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2218 D20 NBG Baseline V4_1_Mod_4_LG STIM.dat';
addpath('z:\sawmillkd\MATLAB\DigitalFilters');
addpath('z:\sawmillkd\MATLAB\DSP');
addpath('z:\sawmillkd\MATLAB\Windows');

[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
BSAG_HEADER = [fp_path '\' fp_file '.txt'];

headerData = ReadBSAGHeader(BSAG_HEADER);

% output file
BSAG_BURSTFILE = [fp_path '\' fp_file '_burst.txt'];
fid_burst = fopen(BSAG_BURSTFILE,'w');

%%
fid = fopen(BSAG_FILE,'rb');

dt_data  = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);
maxEnvelope = sqrt(headerData.LOADING_MAX_VALUE.^2 + headerData.LOADING_MIN_VALUE.^2);
skipval = 0;

nSamples = 2000000;

firstSample = 1   ;

fseek(fid,(firstSample-1)*headerData.DATA_BYTES_PER_UNIT,'bof');


thresholds.hardzero = 100000;
thresholds.pulserise = 6e6;
thresholds.maxysize  = 216477760;
triggering.presamples  = 100;
triggering.postsamples = 5.4e-3/dt_data;



fir.burst=firFilter(0,20000,1/dt_data,256*16,'blackman(x)');

figure(1);
h.burstplot_data = plot(0,0);
h.burstplot_axes = gca;
set(h.burstplot_axes,'YLim',[0 thresholds.maxysize]);
tic;
for k = 1 : floor(headerData.SIGNAL_NUM_SAMPLES/nSamples)
    lastSample = firstSample-1+nSamples;
    disp(lastSample);
    %disp([k firstSample lastSample 100*firstSample / headerData.SIGNAL_NUM_SAMPLES]);
    % read counts as signed int16
    x   = fread(fid, 2*nSamples, '*int16', skipval, 'ieee-be');
    % compute iq magnitude squared
    xm2 = int32(x(1:2:end)).^2+int32(x(2:2:end)).^2;

    % Triggers for bursts
    triggers = TriggerPoints(fir.burst,thresholds.pulserise,triggering.presamples,triggering.postsamples,single(xm2'),0);
     
    
    % Read pulsewidth & prf data for this burst
    for itrig = 1 : length(triggers.lead)
        [pw,prf,prf2,pwr] = GetWaveformProps(thresholds.pulserise,thresholds.hardzero,xm2(triggers.lead(itrig):triggers.trail(itrig)),dt_data);
        ydata = xm2(triggers.lead(itrig):triggers.trail(itrig));
        set(h.burstplot_data,'XData',(0:(length(ydata)-1)),'YData',ydata);
        %drawnow;
        fprintf(fid_burst,'%d %d %d %d %d %d\n',[triggers.lead(itrig)+firstSample-1 triggers.trail(itrig)+firstSample-1 floor(1e9*pw) floor(prf) floor(prf2) floor(pwr)]);
        fprintf(         'PCT=%3d Start: %12d %8dns %6dHz %6dHz Amp: %6d\n',[floor(100*firstSample/headerData.SIGNAL_NUM_SAMPLES) firstSample floor(1e9*pw) floor(prf) floor(prf2) floor(pwr)]);
    end
    
    firstSample=lastSample+1;
end
elapsedTime = toc;
sampPerSec = firstSample / elapsedTime;
disp(['Processing Rate: ' num2str(sampPerSec) ' samples / sec.']);
fclose(fid);
fclose(fid_burst);