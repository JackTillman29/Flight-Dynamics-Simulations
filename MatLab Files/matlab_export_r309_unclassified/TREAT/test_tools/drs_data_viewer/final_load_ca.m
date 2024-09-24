close all;
clear;
clc;

%[a,b]=uigetfile('\\xxx.xxx.xxx.xxx\<name>\EZJ\*.tar');
%[a,b]=uigetfile('D:\EZJA_Files\*.tar');
%total_path=[b a];

%[pathstr,filename,fileext]=fileparts(total_path);
%untar(total_path);
filename = 'D:\EZJ\Coarse_Fine_Analysis\<>';
% in matlab

% debug


dir_list_fine   = dir([filename '\*fine*']);
dir_list_coarse = dir([filename '\*coarse*']);

nfilesc = length(dir_list_coarse);
nfilesf = length(dir_list_fine);
data = [];
dataf = [];

for k = 1 : nfilesc
    sfile = dir_list_coarse(k).name;
    %vals=sscanf(sfile,"channel_and_coarse_analysis_dump_H%d_T%d.txt");
    vals = [0 0]';
    disp(['Processing file: ' sfile]);
    hval=vals(1);
    tval=vals(2);
    temp = load([filename '\' sfile]);
    data = [data; [hval*ones(1024,1) tval*ones(1024,1) temp]];
end

for k = 1 : nfilesf
    sfile = dir_list_fine(k).name;
    %vals=sscanf(sfile,"channel_and_coarse_analysis_dump_H%d_T%d.txt");
    vals = [0 0]';
    disp(['Processing fine file: ' sfile]);
    hval=vals(1);
    tval=vals(2);
    temp = load([filename '\' sfile]);
    dataf = [dataf; [hval*ones(1024,1) tval*ones(1024,1) temp]];
end

idx = find(dataf(:,4) == 0);
dataf(idx,5) = nan;



%% plotting
refIF = 0;
figure;
hold on;
disp('Unique Head Settings:');
disp(unique(data(:,1)));
legstr = {};
disp('Unique Tail Settings:');
disp(unique(data(:,2)));

if((length(unique(data(:,1))) + length(unique(data(:,2)))) == 2)
    singleMode = 1;
else
    singleMode = 0;
end

if(singleMode == 0)
    headFilter = 20;
    tailFilter = [];
else
    headFilter = [];
    tailFilter = [];
end


for kfile = 1 : nfilesc
    iStart = (kfile-1)*1024+1;
    iStop  = iStart + 1024 - 1;
    if(~isempty(headFilter))
        if(data(iStart,1) ~= headFilter)
            continue;
        end
    end
    if(~isempty(tailFilter))
        if(data(iStart,2) ~= tailFilter)
            continue;
        end
    end
    if(singleMode == 0)
        plot(1e-6*(data(iStart:iStop,5)-0*refIF),'.-');
        legstr = {legstr{:},sprintf('H%d\\_T%d',data(iStart,1),data(iStart,2))};
    end
end

if(singleMode == 1)
    %refIF = 1e-3*mean(data(:,5));
    refIF = 670e6*1e-3;
    plot(1e-3*data(:,5)-refIF,'.-');
    title(sprintf('CA Results (rel IF %f)',refIF*1e-3));
    ylabel('kHz');
    
    figure;
    plot(1e-3*dataf(:,5)-refIF,'.');
    title(sprintf('FA Results (rel IF %f)',refIF*1e-3));
    ylabel('kHz');
    
else
    title(['CA Results (IF)']);
    legend(legstr);
ylabel('MHz');
end



grid on;

if(singleMode)
    figure;
    hist(1e-3*data(:,5)-refIF,100);
    xlabel('kHz');
    title(sprintf('Head: %d, Tail %d',hval,tval));   
end