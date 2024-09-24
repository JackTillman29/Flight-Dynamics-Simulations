close all; clear all; clc

N = 101;
N = [101 201];



figure;
for angle_deg = 0:360
    
%     % usage 1
%     N = 101;
%     imagesc(circ_kernel_offset(N,angle_deg))
    
%     % usage 2
%     N = [101 201];
%     imagesc(circ_kernel_offset(N,angle_deg))
    
%     % usage 3
%     M = 100;
%     N = 201;
%     imagesc(circ_kernel_offset(M,N,angle_deg))
    
    % usage 4
    M = 101;
    N = 201;
    dist_from_center = 3;
    circ_rad  = 0.5;
    imagesc(circ_kernel_offset(M,N,angle_deg,dist_from_center, circ_rad))
    
    axis equal
    drawnow
end




