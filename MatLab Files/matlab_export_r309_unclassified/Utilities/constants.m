%% Constants & Converting Functions
k = 1.3805E-23;
c = 2.99792458E8;
r2d = @(x) x.*(180/pi);
d2r = @(x) x.*(pi/180);
m2km = @(x) x.*1e-3;
km2m = @(x) x.*1e3;
hz2khz = @(x) x.*1e-3;
khz2hz = @(x) x.*1e3;
hz2mhz = @(x) x.*1e-6;
mhz2hz = @(x) x.*1e6;
W2dB = @(x)10.*log10(abs(x));
dB2W = @(x) 10.^(x./10);
V2dB = @(x)20.*log10(abs(x));
dB2V = @(x) 10.^(x./20);