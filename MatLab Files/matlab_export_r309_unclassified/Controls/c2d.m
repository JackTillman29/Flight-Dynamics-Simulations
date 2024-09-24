function [numd,dend,k] = c2d(numc,denc,Ts,method,Fwarp)
% This function will convert an s-domain transfer function
% into a z-domain transfer function using the bilinear transform
% It does not require the MATLAB controls toolbox
% K. Sawmiller, 2011
%
% Syntax: [numd,dend,k(optional)] = c2d(numc,denc,Ts,method,option)
% 
% 
% 
% Example: numc = [1.46 1]         (1.46s + 1) / (6.75s + 1)
%          denc = [6.75 1]
%          Ts = 0.01

if(~exist('Fwarp','var') && ~strcmpi(method,'prewarp')) 
    Fwarp = 1.0;
end
if(~exist('method','var'))
    disp('Method: tustin');
    method='tustin';
end

ki = find(numc ~= 0, 1, 'first');
numc = numc(ki:end);
ki = find(denc ~= 0, 1, 'first');
denc = denc(ki:end);

% First, convert to pole-zero form
[zc,pc,k] = zpk(numc,denc);

% Next, convert the s-domain poles/zeros to z-domain
zd = s2z(zc,Ts,method,Fwarp);
pd = s2z(pc,Ts,method,Fwarp);

% Reduced order (n-m) numerator creates (n-m) zeros at infinity
m = length(numc);
n = length(denc);
mismatch = n - m;
if(mismatch > 0)
disp(['Adding ' num2str(mismatch) ' zeros at infinity to account for mismatch']);
end
for k = 1:mismatch
    zd = [zd; -1];
end

% Poles & Zeros are correct. Now convert back to polynomial form.
numd =       poly(zd);
dend =       poly(pd);

% Ensure DC gain is honored
if( (length(denc)-length(numc)) ~= 0)
%dcgc = polyval(numc,0) ./ polyval(denc,0);
%testval =  s2z(1e-3,Ts,'tustin');
%dcgd = polyval(numd,testval) ./ polyval(dend,testval);
%dgain = dcgc / dcgd;
%numd = dgain * numd;
numd = numd * (Ts/2).^(length(denc)-length(numc));
end
numd = numd * numc(1)./denc(1);

% Now convert the z-domain zpk form back to polynomial form
% Apply gain term derived from zpk call
if(nargout == 2)
    % good to go
elseif(nargout == 3)
    k = numd(1);
    numd = numd./k;
else
    error('Incorrect number of output arguments!');
end