close all;
clear all;
clc;

addpath('I:\MATLAB\Visualization');
Fs = 5e4;
t = single(0:1/Fs:1e-3);

yfcn   = @(f,p,t) sum(cos(2*pi*f'*t + p'),1);
yfcnc   = @(f,p,t) abs(sum(exp(1i*(2*pi*f'*t + p')),1));


f = [-2000 -1000 0 1000 2000] + 0*1e4;
p = [0 0 0 0 0];
p = [ ...
       -1.5708
    0.5236
   -3.3161
   -2.7925
   -3.6652 ]';
rstep = 10*pi/180;

hf=figure;
k = 1;subplot(6,1,k);hp(k) = plot(yfcn(f(k),p(k),t));ht(k) = title(sprintf('%1.2f',p(k)));
k = 2;subplot(6,1,k);hp(k) = plot(yfcn(f(k),p(k),t));ht(k) = title(sprintf('%1.2f',p(k)));
k = 3;subplot(6,1,k);hp(k) = plot(yfcn(f(k),p(k),t));ht(k) = title(sprintf('%1.2f',p(k)));
k = 4;subplot(6,1,k);hp(k) = plot(yfcn(f(k),p(k),t));ht(k) = title(sprintf('%1.2f',p(k)));
k = 5;subplot(6,1,k);hp(k) = plot(yfcn(f(k),p(k),t));ht(k) = title(sprintf('%1.2f',p(k)));
subplot(6,1,6);
ysum=yfcnc(f,p,t);
hp(6) = plot(ysum);ht(6) = title(sprintf('Peak: %1.2f',max(ysum)));
ylim([0 5.5]);

set(gcf,'Position',[ 1319          62         584        1038]);

set(gcf,'WindowScrollWheelFcn',@scroll_fcn);


%%

if(0)
    %% Auto - loop
    p = [0 0 0 0 0];
    rstep = 10*pi/180;
    nstep = round((60*pi/180)/rstep);  % appears optimization space is within 60
    % degrees. The pattern repeats using a coarse grid
    
    peakvals = zeros(nstep.^4,5);
    k = 1;
    for k2 = 1 : nstep
        p(2) = (k2-1)*rstep;
        fprintf('%d of %d complete\n',k2,nstep);
        for k3 = 1 : nstep
            p(3) = (k3-1)*rstep;
            for k4 = 1 : nstep
                p(4) = (k4-1)*rstep;
                for k5 = 1 : nstep
                    p(5) = (k5-1)*rstep;
                    
                    peakvals(k,:) = [p(2:5) max(yfcn(f,p,t))];
                    k = k + 1;
                end
            end
        end
    end
    
end