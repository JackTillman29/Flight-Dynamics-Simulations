function png_to_avi2(prefix_name,sort_mode,avi_file,video_profile,video_quality,fps,xyScale)
%filename must be of form <name>_###.png
file_list = dir([prefix_name '_*.png']);
nFiles = length(file_list);

file_num = zeros(1,length(file_list));
for kfile = 1 : nFiles
    this_file_name = file_list(kfile).name;
    file_num(kfile) = sscanf(this_file_name,[prefix_name '_%d.png']);
end

file_num = zeros(1,length(file_list));
for kfile = 1 : nFiles
    this_file_name = file_list(kfile).name;
    file_num(kfile) = sscanf(this_file_name,[prefix_name '_%d.png']);
end

if(strcmp(sort_mode,'ascend'))
    [~,sort_idx]=sort(file_num,'ascend');
else
    [~,sort_idx]=sort(file_num,'descend');
end


figure;

for kfile = 1 : nFiles
    fileToDisplay = file_list(sort_idx(kfile)).name;
    disp( fileToDisplay );
    img = imread(fileToDisplay);
    
    % truncate image
    img = img(:,170:end-170,:); % TODO: generalize to an input
    
    if(xyScale ~= 1)
        img = imresize(img,xyScale);
    end
    if(kfile == 1)
        
        %%
        hi = image(img);
        set(gca,'Visible','off','Position',[0 0 1 1]);
        imX = size(img,2);
        imY = size(img,1);
        set(gcf,'Position',[100 100 imX imY]);
        vw = VideoWriter(avi_file,video_profile);
        if(~strcmp(video_profile,'Uncompressed AVI'))
            vw.Quality = video_quality;
        end
        vw.FrameRate = fps;
        vw.open();
    else
        set(hi,'CData',img);
    end
    drawnow;
    vw.writeVideo(img);
end
close(vw);
    

