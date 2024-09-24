close all;
clear;
clc;
r2d = 180/pi;
d2r = pi/180;

B_a = eye(3)';
u_a = B_a * [1 0 0]';
v_a = B_a * [0 1 0]';
w_a = B_a * [0 0 1]';

yaw = 15*d2r;
pitch = 30*d2r;
roll = 60*d2r;

B_y = RotationMatrixYPR(yaw,0,0)';
u_y = B_y * u_a;
v_y = B_y * v_a;
w_y = B_y * w_a;



B_p = RotationMatrixYPR(0,pitch,0)';
u_p = B_p * u_y;
v_p = B_p * v_y;
w_p = B_p * w_y;

B_r = RotationMatrixYPR(0,0,roll);
u_r = B_r * u_p;
v_r = B_r * v_p;
w_r = B_r * w_p;
