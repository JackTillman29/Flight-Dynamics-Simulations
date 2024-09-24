close all;
clear all;
clc;

% Load file
BSAG_FILE = 'E:\Data from 20150622\BSAG Stim\T2188 S20 NBG Baseline V4.1 Stim.dat';
[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
BSAG_HEADER = [fp_path '\' fp_file '.txt'];
headerData = ReadBSAGHeader(BSAG_HEADER);
fid = fopen(BSAG_FILE,'rb');


bytesPerIQsample = headerData.DATA_BYTES_PER_UNIT;
fseek(fid,0,'eof');
nBytes = ftell(fid);
frewind(fid);
if(nBytes ~= headerData.SIGNAL_NUM_BYTES)
    error('Header bytes does not match file size');
end
return;

%% Purpose of this file is to generate avi files for BSAG data to see what 
% is in the data

fps      = 30;
dt_frame = 1 / fps;
dt_data  = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);

nSamplesPerFrame = floor(dt_frame / dt_data);





sampleHop = 1; % every 100th sample

nIQSamples = nBytes / bytesPerIQsample;
nIsamples = nIQSamples / 2;
nQsamples = nIQSamples / 2;
skipval = (2*sampleHop - 1) * bytesPerIQsample/2;

nFrames = floor(nIQSamples/nSamples);

for k = 1 : nFrames
    disp(k);
    fseek(fid,(k-1)*nSamples*bytesPerIQsample,'bof');
    xi = fread(fid, nSamples, '*int16', skipval, 'ieee-be');
    frewind(fid);
    fseek(fid,(k-1)*nSamples*bytesPerIQsample+2,'bof');
    xq = fread(fid, nSamples, '*int16', skipval, 'ieee-be');

    figure(1);
    plot(abs(xi+0i*xq));
    set(1,'Position',[13 700        1240         218]);
    set(1,'PaperPositionMode','auto');
    print('-dpng','-f1','-r180',['Frame_' num2str(k) '.png']);
end
fclose(fid);