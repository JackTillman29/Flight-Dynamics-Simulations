close all;
clear;
clc;
[a,b]=uigetfile('D:\EZJA_Files\*.tar');
total_path=[b a];

[pathstr,filename,fileext]=fileparts(total_path);
untar(total_path);

% in matlab

% debug


dir_list_fine   = dir([filename '\*fine*']);
%dir_list_coarse = dir([filename '\*coarse*']);

%nfilesc = length(dir_list_coarse);
nfilesf = length(dir_list_fine);
data = [];
for k = 1 : nfilesf
    sfile = dir_list_fine(k).name;
    vals=sscanf(sfile,"channel_and_fine_analysis_dump_H%d_T%d.txt");
    disp(['Processing file: ' sfile]);
    hval=vals(1);
    tval=vals(2);
    temp = load([filename '\' sfile]);
    data = [data; [hval*ones(1024,1) tval*ones(1024,1) temp]];
end

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


for kfile = 1 : nfilesf
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
        plot(1e-6*(data(iStart:iStop,5)-0*refIF),'.');
        legstr = {legstr{:},sprintf('H%d\\_T%d',data(iStart,1),data(iStart,2))};
    end
end

if(singleMode == 1)
    %refIF = 1e-3*mean(data(data(:,5)~=0,5));
    plot(data(:,5),'.');
    title(sprintf('FA Results'));
    ylabel('Hz');
    xlimits = xlim();
    h1 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[500 550 550 500]*1e6,'r');
    h2 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[550 600 600 550]*1e6,'r');
    h3 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[600 650 650 600]*1e6,'r');
    h4 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[650 700 700 650]*1e6,'r');
    h5 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[700 750 750 700]*1e6,'r');
    h6 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[750 800 800 750]*1e6,'r');
    set([h1 h2 h3 h4 h5 h6],'EdgeColor','k','FaceColor','k','FaceAlpha',0.1);
else
    title(['FA Results (IF)']);
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