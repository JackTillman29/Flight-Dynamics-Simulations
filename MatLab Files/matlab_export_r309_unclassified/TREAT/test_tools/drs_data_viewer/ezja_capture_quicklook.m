close all;
%clear;
clc;
% TODO
% 1 - Pulse present threshold is hard coded in cb_iq.m
% 2 - Clean up Hilbert conversion
%

%drs_decimation = 16;
%Fs = 100e6/drs_decimation;
%Ts = single(1/Fs);

matlab_root = 'I:\MATLAB';
addpath(matlab_root);
addpath([matlab_root '\Analysis']);
addpath([matlab_root '\Printing']);
addpath([matlab_root '\Utilities']);
pp=PrepForPrint;

ezja_header = 1;


% Function to convert sample time (s) to byte in file
% Arguments: time (s), Fs (Hz), Bits/Sample (#)
sampleTimeToByte = @(t,Fs,bps) floor(t * Fs * bps/8);

if(~exist('reload_flag'))
    reload_flag = 0;
end

if(~exist('lastfile'))
    if(ezja_header == 0)
        lastfile = '*.tmp';
    else
        lastfile = '*.dat';
    end
end

%if(reload_flag == 0)

% === Setup Defaults ===
if(~exist('ref_channel'))
    ref_channel = 1;
end

% if(~exist('drs_decimation'))
%     drs_decimation = 16;
% end

if(~exist('ref_channel_dop_offset'))
    ref_channel_dop_offset = 0; %17379-500;
end

if( ~exist('lastpath') )
    lastpath = '<enter default path>';
end
% === End Setup Defaults ===

[a,b]=uigetfile([lastpath lastfile],'Select Channel 1');
inChanStr = [b a];
lastpath = b;
lastfile = a;



%[a,b]=uigetfile([lastpath '*.tmp'],'Select Channel 1');
% Replace _1.tmp in first file choice with _2.tmp
a=strrep(a,'_1.tmp','_2.tmp');
outChanStr = [b a];

    


voltsPerBit = 1.0;
%drs_decimation = 16;




if(ezja_header == 0)
    header1 = parse_midas_header(inChanStr);
    header2 = parse_midas_header(outChanStr);
else
    header1.data_start = 1;
    header2.data_start = 1;
end
fid1 = fopen(inChanStr,'rb');
fid2 = fopen(outChanStr,'rb');

if(ezja_header == 0)
    % Get sample time & rate
    fseek(fid1,264,'bof');
    Ts = fread(fid1,1,'real*8');
    Fs = 1.0 / Ts;
    fseek(fid1,0,'bof');
else
    Fs = 2.6e9;
    Ts = 1.0 / Fs;
    bps = 12;
end


% Just get the total number of bytes in the file
fseek(fid1,0,'eof');        % move to end of file
nTotalBytes = ftell(fid1);   % capture file pointer (number of bytes)
fseek(fid1,0,'bof');        % move to beginning of file



if(ezja_header == 0)
    nSamplesInFile = (nTotalBytes - 512) / 2.0 / 2.0; % two bytes per sample, I/Q samples
else
    nSamplesInFile = floor(nTotalBytes * 8 / 12);
end


if(~exist('timeArray'))
    timeArray = [0 nSamplesInFile]./Fs; % Amount of time (sec) in file
end

fseek(fid1,header1.data_start-1,'bof');
fseek(fid2,header2.data_start-1,'bof');


res = inputdlg({ ...
    sprintf( ...
    'Time Window [%5.1f - %5.1f]:',0,nSamplesInFile./Fs), ...
    'Ref Channel:', ...
    'Dop Offset (Hz):', ...
    'Preview Decimation:'},'?',1, ...
    {sprintf('%8.4f,%8.4f',timeArray(1),timeArray(2)),num2str(ref_channel),num2str(ref_channel_dop_offset),'50'});

prev_dec = str2num(res{4});
% Check for brackets
temp = res{1};
temp = erase(temp,'[');
temp = erase(temp,']');
timeWindow = sscanf(temp,'%f,%f')';
% end check for brackets

ref_channel = sscanf(res{2},'%d');


ref_channel_dop_offset = sscanf(res{3},'%f')';
%drs_decimation = sscanf(res{4},'%f');

%else
%    prev_dec = 1;
%    timeWindow = reload_window;
%end
%byteStart = Fs * timeWindow(1) * 2 * 2;
byteStart = sampleTimeToByte(timeWindow(1),Fs,bps);
disp(['Byte Start: ' num2str(byteStart)]);
if(ezja_header == 0)
    Nsamples = diff(timeWindow) * Fs * 2;
else
    Nsamples = diff(timeWindow) * Fs;
end
Nsamples = round(Nsamples);
if(mod(Nsamples,2))
    Nsamples = Nsamples + 1;
end
fseek(fid1,byteStart,'cof');
fseek(fid2,byteStart,'cof');
disp(['Loading DRS Data ' inChanStr]);
disp(['Loading DRS Data ' outChanStr]);
if(ezja_header == 0)
    data1 = fread(fid1,Nsamples,'int16=>int16','ieee-le');
    data2 = fread(fid2,Nsamples,'int16=>int16','ieee-le');
else
    data1 = fread(fid1,Nsamples,'bit12=>int16','ieee-le');
    data2 = fread(fid2,Nsamples,'bit12=>int16','ieee-le');
end

if(ezja_header == 0)
    dat1.pageMemory = [data1(1:2:end) data1(2:2:end)];
    dat2.pageMemory = [data2(1:2:end) data2(2:2:end)];
else
    dat1.pageMemory = [data1];
    dat2.pageMemory = [data2];
end

disp(['Byte Stop: ' num2str(ftell(fid1))]);
fclose(fid1);
fclose(fid2);
disp('Done');




xvec = 1e6*Ts * (0:length(dat1.pageMemory(:,1))-1);
h.figure = figure;

h.a1 = subplot(2,1,1);
h.a2 = subplot(2,1,2);
plot_str = '-';
if(ezja_header == 0)
h.line1 = plot(h.a1, ...
    xvec(1:prev_dec:end), ...
    10*log10((single(dat1.pageMemory(1:prev_dec:end,1)).^2+single(dat1.pageMemory(1:prev_dec:end,2)).^2)) ...
    ,plot_str);
else
    dat1.temp = single(dat1.pageMemory(1:prev_dec:end));
    dat2.temp = single(dat1.pageMemory(1:prev_dec:end));
    dat1.tempf = fft(dat1.temp);
    dat2.tempf = fft(dat2.temp);
    dat1.tempf(floor(length(dat1.tempf)/2):end) = 0;
    dat2.tempf(floor(length(dat2.tempf)/2):end) = 0;
    dat1.tempc = ifft(dat1.tempf);
    dat2.tempc = ifft(dat2.tempf);
    
    h.line1 = plot(h.a1, ...
    xvec(1:prev_dec:end), ...
    20*log10(abs(dat1.tempc)) ...
    ,plot_str);
end
%set(gca,'UserData',1:length(dat1.pageMemory(:,1)));

if(ezja_header == 0)
h.line2 = plot(h.a2, ...
    xvec(1:prev_dec:end), ...
    10*log10((single(dat2.pageMemory(1:prev_dec:end,1)).^2+single(dat2.pageMemory(1:prev_dec:end,2)).^2)) ...
    ,plot_str);
else
    h.line1 = plot(h.a2, ...
    xvec(1:prev_dec:end), ...
    20*log10(abs(dat2.tempc)) ...
    ,plot_str);
end
%set(gca,'UserData',1:length(dat1.pageMemory(:,1)));

subplot(2,1,1);
xlabel('Time (\mus)');
ylabel('20log_{10} Counts');
title(sprintf('Channel 1 [%5.2f to %5.2f]',timeWindow(1),timeWindow(2)));

subplot(2,1,2);
xlabel('Time (\mus)');
ylabel('20log_{10} Counts');
title(sprintf('Channel 2 [%5.2f to %5.2f]',timeWindow(1),timeWindow(2)));
linkaxes([h.a1 h.a2],'x');

hfig = gcf;
hcmenu = get(hfig,'uicontextmenu');
% figure doesn't contain a context menu object, create a new one.
if(isempty(hcmenu))
    hcmenu = uicontextmenu;
    item1 = uimenu(hcmenu, 'Label', 'Show Spectrum', 'Callback', 'cb_spectrum');
    item2 = uimenu(hcmenu, 'Label', 'Show RDM', 'Callback', 'cb_rdm');
    item3 = uimenu(hcmenu, 'Label', 'Show I/Q', 'Callback', 'cb_iq');
    item4 = uimenu(hcmenu, 'Label', 'Set Ref CH1', 'Callback', 'cb_set_ref1');
    item5 = uimenu(hcmenu, 'Label', 'Set Ref CH2', 'Callback', 'cb_set_ref2');
    item6 = uimenu(hcmenu, 'Label', 'Zoom/Resample', 'Callback', 'cb_get_time');
    item7 = uimenu(hcmenu, 'Label', '** Make Bad News Plot **', 'Callback', 'cb_plot_stft');
    item8 = uimenu(hcmenu, 'Label','Frequency Mod','Callback','cb_fmod');
    set(hfig,'UIContextMenu',hcmenu);
end

if(ref_channel == 1)
    cb_set_ref1();
else
    cb_set_ref2();
end

pauseThreshold = 2e7;
add_analysis_callbacks;

set(gcf,'KeyPressFcn',@drs_key_press);