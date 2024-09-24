% close all; clear all; clc

% nt = 4;
% f0 = 10e3;
% df = 500;

fv = [0:df:(nt-1)*df] + f0;



Fs = 30*max(fv);
dt = 1/Fs;

t = [0:dt:(2/min(fv) - dt)];

siggen = @(f,p,t) sum(cos(2*pi*f'*t + p'),1);

% nrep = 1e6;

peakval = zeros(1,nrep);
meanval = zeros(1,nrep);
for irep = 1:nrep
    
    if(mod(irep,floor(nrep*0.10)) == 0)
        clc
        % disp(['irep = ',num2str(irep),' / ',num2str(nrep)])
        disp([num2str(irep/nrep * 100),' % done'])
        
%         if(~exist('hfig'))
%             hfig = figure;
%             hp = plot(x);
%         else
%             set(hp,'YData',x)
%         end
%         drawnow;
    end
    
    pv = [0 2*pi*rand(1,nt-1)];
    x = siggen(fv,pv,t) ./ length(fv);
    
    peakval(irep) = max(abs(x));
    meanval(irep) = mean(abs(x));
    
end

mpr = meanval ./ peakval;

% nbins = 10000;
% [hmpr,empr] = histcounts(mpr,nbins);
% cmpr = empr(2:end) - diff(empr);

% figure;
% subplot(2,1,1)
% plot(mpr,'.')
% 
% subplot(2,1,2)
% plot(cmpr,hmpr)
% 


