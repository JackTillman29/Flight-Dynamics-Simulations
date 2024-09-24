function png_to_agif(prefix_name,sort_mode,gif_file,delay_s,xyScale)
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
    else
        set(hi,'CData',img);
    end
    drawnow;
    
    [A,map] = rgb2ind(img,256);
    if(kfile == 1)
        imwrite(A,map,gif_file,'gif','LoopCount',Inf,'DelayTime',delay_s);
    else
        imwrite(A,map,gif_file,'gif','WriteMode','append','DelayTime',delay_s);
    end
    
end

    

