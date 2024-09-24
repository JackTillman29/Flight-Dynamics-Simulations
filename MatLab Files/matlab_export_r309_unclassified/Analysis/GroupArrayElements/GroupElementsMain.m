close all; clear all; clc


figure
add_group_tool

N = 100
scale = 5;
x = scale*(2*rand(1,N)-1);
y = scale*(2*rand(1,N)-1);
plot(x,y,'b.')


