close all; clear all; clc

large_mat = zeros(100,100);

small_mat = ones(3,4);
small_mat = rand(70) > 0.3;


figure;
istep = 0;
radius = linspace(1,3,100);
for th = linspace(0,2*pi,100)

    istep = istep + 1;
    idx = 0.5*radius(istep)*size(large_mat,1) + size(large_mat,1)/(2*radius(istep)*cos(th));
    idy = 0.5*radius(istep)*size(large_mat,2) + size(large_mat,2)/(2*radius(istep)*sin(th));
%     idx = idx .* radius(istep);
%     idy = idy .* radius(istep);
    new_mat = embedMatIntoMat(small_mat, large_mat, round([idx idy]));
    imagesc(new_mat.')
    axis equal;
    hold on;
    plot(idx,idy,'go','MarkerFaceColor','g','MarkerSize',10)
    hold off;
    drawnow;
end


