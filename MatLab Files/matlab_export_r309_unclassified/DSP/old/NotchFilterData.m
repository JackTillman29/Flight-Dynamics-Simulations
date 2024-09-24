function DataOut = NotchFilterData(Data,Ts,Fn,AttenDb,Zeta)
% Notch filter design code credited to  Tony Nguyen.
% Notch filter is of the form:
% K*(1 + a1*z^-1 + a2*z^-2)/(1 + b1*z^-1 + b2*z^-2)

%%% Inputs:
%Ts = 0.001;         % sampling time in seconds
%f  = 37.0			% frequency in Hz
%db = 20;            % magnitude in dB
%zN = 0.01;          % damping ratio in the numerator


f = Fn;
wd = 2*pi*f;

db = AttenDb;
zN = Zeta;

wa = 2/Ts * tan(wd*Ts/2);
f = wa/(2*pi);

temp = 10^(db/20);
zD = zN*temp;
w = 2*pi*f;
%num = [1 2*zN*w w*w];
%den = [1 2*zD*w w*w];

% Discrete notch filter using bilinear transform
% s = 2/T*(z-1)/(z+1)
n1 =  4 + 4*zN*w*Ts + w*w*Ts*Ts;
n2 = -8 + 2*w*w*Ts*Ts;
n3 =  4 - 4*zN*w*Ts + w*w*Ts*Ts;
d1 =  4 + 4*zD*w*Ts + w*w*Ts*Ts;
d2 = -8 + 2*w*w*Ts*Ts;
d3 =  4 - 4*zD*w*Ts + w*w*Ts*Ts;
numd = [n1 n2 n3];
dend = [d1 d2 d3];
0.25*[numd ; dend];	% divide both num and den by 4
% Put in the form:
%    K*(1 + a1*z^-1 + a2*z^-2)/(1 + b1*z^-1 + b2*z^-2)
kf  = n1/d1;
a1 = n2/n1;
a2 = n3/n1;
b1 = d2/d1;
b2 = d3/d1;

    for k = 3:length(Data)
        Data(k) = -b1*Data(k-1) + -b2*Data(k-2) + kf*(Data(k) + a1*Data(k-1) + a2*Data(k-2));
    end
    DataOut = Data;
end