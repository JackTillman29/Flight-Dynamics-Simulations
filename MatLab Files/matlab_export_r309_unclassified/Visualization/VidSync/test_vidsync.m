close all;
clear all;
clear classes;

mainFig = figure;
masterFps = 30;

vid.compressionMethod = 'H.264';% 'None', 'H.264', 'Motion JPEG', or 'Motion JPEG 2000'.

a = vidSyncObj('fort.2000.avi',mainFig,masterFps);
b = vidSyncObj('fort.2002.avi',mainFig,masterFps);
c = vidSyncObj('fort.2004.avi',mainFig,masterFps);

% note: frame rate factors must be integers at this time
% a value of 1 means draw every frame at master frame rate
% a value of 2 means draw every frame at half (1/2) of master 
% a value of 5 = 1/5 of master frame rate

a.setup( ...
    [0 0.1 0.7 0.7], ...  % start x,y, fraction x,y on figure
    1);                      % frame rate factor
b.setup([0.5 0.0  0.5 0.5],1);
c.setup([0.5 0.5  0.5 0.5],1);

set(mainFig,'Color','k');
set(mainFig,'Position',[0 0 1280 720]);

nFrames = min([a.nMasterFrames b.nMasterFrames c.nMasterFrames]);

outputVidObj = VideoWriter('output','MPEG-4');
outputVidObj.Quality = 95;

open(outputVidObj);
disp('Draw window to appropriate size, then press space bar');
pause;
for k = 1 : nFrames
    outputVidObj.writeVideo(getframe(mainFig));
    a.draw();
    b.draw();
    c.draw();
end

close(outputVidObj);


