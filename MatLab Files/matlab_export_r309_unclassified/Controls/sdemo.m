% z transform

% create s grid
[res,ims]=meshgrid(-10:10,-10:10);

SP=res+1i*ims;

Ts = 0.001;
ZP = (2+SP.*Ts) ./ (2-SP.*Ts);