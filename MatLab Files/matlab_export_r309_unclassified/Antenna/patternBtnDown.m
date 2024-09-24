cp = get(gca,'CurrentPoint');
azSteer = cp(1,1)*pi/180;
elSteer = cp(1,2)*pi/180;

[uxs,uys,uzs]=sph2cart(azSteer,elSteer,1);
phs_steer = -2*pi / lambda * p * [uxs uys uzs]';
phsnz = phsnz_mag*pi/180*randn(Ne,1);

ws = exp(1i*(phs_steer + phsnz));

w = wa .* ws;
PW = w.' * exp(1i*P) ./ sqrt(Ne);
Res = reshape(PW,nEl,nAz);

set(hi,'CData',20*log10(abs(Res)));
sEarth.CData = [get(hi,'CData') 0*get(hi,'CData')-99];