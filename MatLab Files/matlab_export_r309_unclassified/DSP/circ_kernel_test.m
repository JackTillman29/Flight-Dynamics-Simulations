close all; clear all; clc

figure;
for setN = 5:400
    N = [setN];      % circular
    N = [setN setN]; % circular
%     N = [100 setN];  % elliptical
    
    imagesc(circ_kernel(N))
    axis equal
    drawnow;
end

