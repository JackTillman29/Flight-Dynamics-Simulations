close all; clear all; clc

N = [100 100];
x = rand(N);

ix = 1:N(1);
iy = 1:N(2);
[ixg,iyg] = meshgrid(ix,iy);

kernel = ones(5);
kernel = kernel ./ sum(kernel(:));
PWR = conv2(x,kernel,'same');

det_thresh = 0.6;
A = PWR > det_thresh;
centroids = detectionCentroid2d(A.',PWR.');

figure;
subplot(1,2,1);
imagesc(PWR.');
hold on;
plot(centroids.centroid_ind(:,1),centroids.centroid_ind(:,2),'k.')

subplot(1,2,2);
imagesc(A.');
hold on;
plot(centroids.centroid_ind(:,1),centroids.centroid_ind(:,2),'r.','MarkerSize',10)
colormap(gca,gray)




