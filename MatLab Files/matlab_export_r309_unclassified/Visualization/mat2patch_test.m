close all; clear all; clc



kernel = zeros(5,5);
kernel(3,1) = 1;
kernel(2:4,2) = 1;
kernel(1:5,3) = 1;
kernel(2:4,4) = 1;
kernel(3,5) = 1;

C = rand(100,200);
C = conv2(C,kernel,'same')./sum(kernel(:));

figure;
imagesc([0 1],[0 1],C); hold on;

imagesc([0 1],[2 3],C); hold on;
ylim([0 4])

% figure; imagesc(C)
% return

alpha_map = ones(size(C));
alpha_map(C < 0.5) = 0;
[px,py,pc,palpha] = mat2patch(C,alpha_map);

C2 = rand(100,200);
C2 = conv2(C2,kernel,'same')./sum(kernel(:));
alpha_map2 = ones(size(C2));
alpha_map2(C2 < 0.5) = 0;
[px2,py2,pc2,palpha2] = mat2patch(C2,alpha_map2);


cm = jet(100);
datlims = [min(pc) max(pc)];
ncmPts = size(cm,1);
pcm_ind = zeros(ncmPts,1);
pcm_ind = round(interp1(datlims,[1 ncmPts],pc,'linear'));
pcm = reshape(cm(pcm_ind,:),[length(pcm_ind),1,3]);

cm = gray(100);
datlims = [min(pc2) max(pc2)];
ncmPts = size(cm,1);
pcm_ind2 = zeros(ncmPts,1);
pcm_ind2 = round(interp1(datlims,[1 ncmPts],pc2,'linear'));
% pcm2 = zeros(ncmPts,1,3);
pcm2 = reshape(cm(pcm_ind2,:),[length(pcm_ind2),1,3]);

figure
hp = patch(px,py,pc);
set(hp,'edgecolor','none','FaceAlpha','flat','FaceVertexAlphaData',palpha)
applyColormapToPatch(jet(100),hp)

hold on;

hp2 = patch(px2+1.5,py2,pc2);
set(hp2,'edgecolor','none','FaceAlpha','flat','FaceVertexAlphaData',palpha2)
applyColormapToPatch(parula(100),hp2)


for k = 1:5
    applyColormapToPatch(gray(100),hp); pause(0.2);
    applyColormapToPatch(cool(100),hp); pause(0.2);
    applyColormapToPatch(hot(100),hp);  pause(0.2);
    
    clims = [0.5 0.9];
    applyColormapToPatch(gray(100),hp,clims); pause(0.2);
    applyColormapToPatch(cool(100),hp,clims); pause(0.2);
    applyColormapToPatch(hot(100),hp,clims);  pause(0.2);
    
    clims = [0.5 0.6];
    applyColormapToPatch(gray(100),hp,clims); pause(0.2);
    applyColormapToPatch(cool(100),hp,clims); pause(0.2);
    applyColormapToPatch(hot(100),hp,clims);  pause(0.2);
end










