function png_to_avi(search_name,video_file,video_profile,video_quality)
%search_name should return from "dir(search_name)" the sequence of the
%video frames. So, for example: search_name = 'cp_*_JS_010dB.png'.
% Note: to write out frames with leading zeros, use %04d format (leading 0)
video_quality = video_quality;
file_list = dir(search_name);
nFiles = length(file_list);

file_num = zeros(1,length(file_list));
for kfile = 1 : nFiles
    this_file_name = file_list(kfile).name;
    %file_num(kfile) = sscanf(this_file_name,[prefix_name '_%d.png']);
    file_num(kfile) = kfile;
end

[~,sort_idx]=sort(file_num);

vw = VideoWriter(video_file,video_profile);
figure;

for kfile = 1 : nFiles
    fileToDisplay = file_list(sort_idx(kfile)).name;
    disp( fileToDisplay );
    img = imread(fileToDisplay);
    if(kfile == 1)
        hi = image(img);
        set(gca,'Visible','off','Position',[0 0 1 1]);
        imX = size(img,2);
        imY = size(img,1);
        set(gcf,'Position',[100 100 imX imY]);
        vw.Quality = video_quality;
        vw.open();
    else
        set(hi,'CData',img);
    end
    drawnow;
    vw.writeVideo(img);
end
close(vw);
