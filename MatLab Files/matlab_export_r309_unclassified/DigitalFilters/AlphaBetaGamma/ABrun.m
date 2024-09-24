a = 0.555;
Ts = 61.44e-3;
b = 0.022;

d = load('fort.1006');
ttrue = d(:,1);
rtrue = d(:,2);
rinit = d(1,3);
rdinit = d(1,4);

%t = ttrue(1):Ts:ttrue(end);
%dr = (rtrue(end)-rtrue(1))/(ttrue(end)-ttrue(1));
%rtrue = rtrue(1)+t*dr;


err = 0.0*t;

nFrames = length(ttrue);

PHI = [1 Ts;0 1];
K = [a b/Ts]';
xbar = zeros(2,nFrames);
xhat = xbar;
% initialize states
k = 1;
xbar(:,k) = [rinit rdinit]';
err(1) = 0.0;


for k = 2 : nFrames
    DTVAL = ttrue(k) - ttrue(k-1);
    K(2) = b/DTVAL;
    PHI = [1 DTVAL;0 1];
    disp(DTVAL)
    xhat(:,k) = PHI*xbar(:,k-1);
    err(k) = rtrue(k) - xhat(1,k);
    xbar(:,k) = xhat(:,k) + K * err(k);
end