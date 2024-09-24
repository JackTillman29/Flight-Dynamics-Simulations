function f = abgFilter(A,B,G,Ts)
% K. Sawmiller
FLT.Phi = [1 Ts 0.5*Ts*Ts;0 1 Ts;0 0 1];
FLT.xhat = [0 0 0]';
FLT.xbar = [0 0 0]';
FLT.gain = [A B G]';

f = class(FLT,'abgFilter');

end