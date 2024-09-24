close all; clear all; clc


nx = 100;
ny = 200;

x = linspace(0,1,nx);
y = linspace(5,10,ny);
z = zeros(ny,nx);
z(:,1:nx/2) = 1;

z = rand(ny,nx);


figure; add_analysis_callbacks;
imagesc(x,y,z); set(gca,'YDir','normal')
colorbar
title('Try the "Grab Slice" tools on the toolbar')
xlabel('Doppler [bin]')
ylabel('Range [bin]')


% % figure; plot(rand(1,100)); add_analysis_callbacks
% % 
% % figure; imagesc(rand(16))
% % 
% % 
% % imwrite(randi([1,100],16),parula(100),'random_colors.png')
