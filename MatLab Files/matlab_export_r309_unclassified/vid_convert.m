close all;
clear all;
clc;

vid_in = 'Z:\path\to\video\cool_video.avi';
[a,b]=fileparts(vid_in);


Quality = 40;

vid_out100 = [a '\' b '_Q' num2str(Quality) '.avi'];

%VideoWriter.getProfiles

vin = VideoReader(vid_in);
vout100 = VideoWriter(vid_out100,'Motion JPEG AVI');
vout100.Quality = Quality;

timeFactor = [ ...
    0.0 1
    1.0 1
    1.1 2
    999 2];

tmax = 999;


open(vout100);
disp(['Processing...']);

frCounter = 1;

while(1)
    draw = 1;
%     N = interp1(timeFactor(:,1),timeFactor(:,2),vin.CurrentTime,'nearest');
%     if(frCounter == N) 
%         draw = 1;
%         frCounter = 1;
%     else
%         draw = 0;
%         frCounter = frCounter + 1;
%     end
    
    if(draw)
    
    try
        fr = vin.readFrame();
    catch
        disp('no more input frames!!');
        break;
    end
    
    end
    
    vout100.writeVideo(fr);
    if(vin.CurrentTime > tmax)
        close(vout100);
        disp('timed out.');
        return;
    end
end

%close(vout60);
%close(vout80);
%close(vout90);
close(vout100);

