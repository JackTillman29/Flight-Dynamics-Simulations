close all;
clear all;

addpath('C:\Users\<user>\Documents\MATLAB\MATLAB\Printing\');

dt = 0.01;

t = 0:dt:1;
f = 2;

figure;
plot3(t,cos(2*pi*f*t),0*t);
hold on;
plot3(t,0*t,sin(2*pi*f*t));
hold off;

hp = patch( ...
    2*[0 0 0 0], ...
    2*[0 0 1 1]-1, ...
    2*[0 1 1 0]-1, ...
    'r');
set(hp,'EdgeColor',0.4*[1 1 1], ...
    'FaceColor',0.8*[1 1 1], ...
    'FaceAlpha',0.3);

