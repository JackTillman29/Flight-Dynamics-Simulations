close all;
clear;
clc;

ub = [1 0 0]'; % forward
vb = [0 1 0]'; % right
wb = [0 0 1]'; % down

R_I2B = RotationMatrixYPR(15*pi/180,30*pi/180,40*pi/180);
R_B2I = R_I2B';

% expect forward to now be -down
ui = R_B2I * ub;
vi = R_B2I * vb;
wi = R_B2I * wb;

% now reconstruct R matrix using basis set

R_I2B_uv1 = [ ...
    ub'*ui vb'*ui wb'*ui
    ub'*vi vb'*vi wb'*vi
    ub'*wi vb'*wi wb'*wi ];

R_I2B_uv2 = [ub vb wb] * [ui vi wi]';
R_B2I_uv2 = [ub vb wb] * [ui vi wi];