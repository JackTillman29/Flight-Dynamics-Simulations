numc = [52.2  1800];
denc = [0.079 1];
Ts = 0.0625;

zc = roots(numc)
pc = roots(denc)

% pole/zero map
zd = exp(zc*Ts)
pd = exp(pc*Ts)

numd = poly(zd)
dend = poly(pd)

% dc gain of input
dcgc = polyval(numc,0) ./ polyval(denc,0)

dcgd = polyval(numd,1) ./ polyval(dend,1)

dgain = dcgc / dcgd

numd = dgain * numd

dcgd = polyval(numd,1) ./ polyval(dend,1);

% works!

% now check a 1st order num/2nd order den
clear;clc;
numc = [5 250];
denc = [2 5 250];
Ts = 0.0625;

zc = roots(numc)
pc = roots(denc)

% pole/zero map
zd = exp(zc*Ts)
pd = exp(pc*Ts)

% append order mismatch zeros
m = length(numc);
n = length(denc);
mismatch = n - m
for k = 1:mismatch
    zd = [zd -1];
end

numd = poly(zd)
dend = poly(pd)

% dc gain of input
dcgc = polyval(numc,0) ./ polyval(denc,0)

dcgd = polyval(numd,1) ./ polyval(dend,1)

dgain = dcgc / dcgd

numd = dgain * numd

dcgd = polyval(numd,1) ./ polyval(dend,1);

% works!

% try tustin
%%
clear;clc;
numc = [5 250];
denc = [2 5 250];
Ts = 0.0625;

zc = roots(numc)
pc = roots(denc)

% pole/zero map
zd = s2z(zc,Ts,2*pi)
pd = s2z(pc,Ts,2*pi)

% append order mismatch zeros
% m = length(numc);
% n = length(denc);
% mismatch = n - m
% for k = 1:mismatch
%      zd = [zd -1];
% end

numd = poly((zd))
dend = poly((pd))

% dc gain of input
dcgc = polyval(numc,0) ./ polyval(denc,0)

dcgd = polyval(numd,1) ./ polyval(dend,1)

dgain = dcgc / dcgd

numd = dgain * numd

dcgd = polyval(numd,1) ./ polyval(dend,1);

