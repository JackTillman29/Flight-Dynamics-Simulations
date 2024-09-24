close all;
clear all;
clc;
hf2 = figure;
hax1 = axes();
set(hax1,'NextPlot','add');

figure('position',[520 121 996 677]);
haxhp1 = subplot(2,1,1)
hp = plot(0,0);

yfcn   = @(fc,f,p,t) 0.5*cos(2*pi*fc(1)*t) + sum(cos(2*pi*fc'*t).*cos(2*pi*f'*t + p'),1);

fv = [1:1:5]*1000;
%fv = [500 1000];
% fv = [500 1000 1500];
fc = 2*fv(end) * ones(1,length(fv));
pv = 0*fv+2*pi*rand(1,length(fv))
ov = 0*fv;           % changes to "1" if this phase has been involved in optimization

d1v = 0*fv;
d2v = 0*fv;
phase_opt = 180*pi/180;
%
N = length(fv);

ylim([-1 1]*N);
dt = 1e-6; Fs = 1/dt;
t=[0:dt:(10/min(fv))];
hL = line([t(1) t(end)],[0 0],'Color','r');

haxhp2 = subplot(2,1,2);
set(haxhp2,'NextPlot','add')
% hp2 = plot(0,0,'r.-');
hp3 = plot(1e6,1e6,'.');



for kopt = 1 : 1e6
    pv = mod(pv - pv(1),2*pi);
    
    yout = yfcn(fc,fv,pv,t);
    set(hp,'XData',t,'YData',yout);
    yrms = sqrt(mean((yout./max(abs(yout))).^2));
    set(hL,'YData',[1 1]*yrms);
    
    fyout = fftshift(fft(yout)) / length(yout);
    df = Fs/length(yout);
    farr = Fs/2 * linspace(-1,1,length(yout));
    minmaxF = 10e3;
%     fidx = find(abs(farr) < minmaxF);
    
    set(hp3,'color',0.7*ones(1,3),'linewidth',1)
    %hp3 = plot(haxhp2,farr,10*log10(abs(fyout)),'b.-','linewidth',2);
    
%     set( hp2,'XData',farr,'YData',10*log10(abs(fyout)) )
    xlim(haxhp2,[-1 1]*minmaxF)
   % plot(hax1,kopt,pv','b.');
    
    
    
    % find peak
    [p,pj]=max(abs(yout));
    
    % look at impacts for this point
    for k = 1 : length(fv)
        pv2 = pv;
        pv2(k) = pv2(k) + phase_opt;
        %yout2(k) = abs(yfcn(fv,pv2,t(pj))) - p;
        yout2(k) = max(abs(yfcn(fc,fv,pv2,t))) - p;
    end
    
    % find most negative
    [impmin,impminj] = min(yout2);

    
    pv(impminj) = pv(impminj) + phase_opt;
    
   
    ov(impminj) = 1;
    % once only a single phase value remains untouched.. this seems to be
    % the optimal solution. This triggers after the last untouched phase is
    % perterbed.
    title(haxhp1,{num2str(pv,'%5.3f '),num2str(20*log10(1./max(abs(yout))))})
    if(sum(ov) == N)
        disp('Total Loss (dB)');
        disp(20*log10(1./max(yout)))
        %break;
    end
    %pv(d1i) = pv(d1i) + .01;
    
    
    
    drawnow;
    

end

return;

%% Test
figure;
hold on;
atten = 1./max(yout);
sumwav = 0*t;
for k = 1 : N
    thiswav=atten*cos(2*pi*fv(k)*t+pv(k));
    plot(t,thiswav);
    sumwav = sumwav + thiswav;
end

plot(t,sumwav,'Color',0.6*[1 1 1]);

