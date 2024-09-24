close all;
clear;
clc;

% create the impulse response (this is the unknown signal)
u = [1 2 3 4 5];

% create the known signal
s = [1 -1 2];

% create the "as measured" signal
y = conv(u,s,'full');


C = convmatrix(s,u);
disp('Convolution Matrix');
disp(C);

disp('conv() result');
disp(y);

disp('C * u = y result');
disp((C*u')');

disp('Deconvolution - Pseudoinverse');
u1 = (inv(C'*C)*C')*y'

% SVD approach
disp('Deconvolution - SVD/PCA');
InvC = pinv_svd(C);
u2 = InvC * y'