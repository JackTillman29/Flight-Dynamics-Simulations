close all;
clear all;
clc;

figure;
hp = plot(0,0);

yfcn   = @(f,p,t) cos(2*pi*f*t + p);
ydfcn  = @(f,p,t) -sin(2*pi*f*t+p) * (2*pi*f);
yddfcn = @(f,p,t) (2*pi*f) * (-cos(2*pi*f*t+p) * (2*pi*f) );

fv = [-2 -1 0 1 2]*1000+2000;
pv = [0 0 0 0 0]*pi/180;
ov = [0 0 0 0 0];           % changes to "1" if this phase has been involved in optimization

d1v = 0*fv;
d2v = 0*fv;


N = length(fv);

ylim([-1 1]*N);
t=[0:1e-6:(1/1000)];
hL = line([t(1) t(end)],[0 0],'Color','r');
for kopt = 1 : 1e6
    yout = 0;
    for k = 1 : N
        yout = yout + yfcn(fv(k),pv(k),t);
    end
    set(hp,'XData',t,'YData',yout);
    yrms = sqrt(mean((yout./max(abs(yout))).^2));
    set(hL,'YData',[1 1]*yrms);
    
    % find peak
    [p,pj]=max(abs(yout));
    
    % look at derivatives for this point
    for k = 1 : length(fv)
        d1v(k) = ydfcn(fv(k),pv(k),t(pj));
        d2v(k) = yddfcn(fv(k),pv(k),t(pj));
    end
    
    if(yout(pj) > 0) % maximum
        % find derivative which is most NEGATIVE
        [d1,d1i]=min(d1v);
    else % minimum
        % find derivative which is most POSITIVE
        [d1,d1i]=max(d1v);
    end
    ov(d1i) = 1;
    % once only a single phase value remains untouched.. this seems to be
    % the optimal solution. This triggers after the last untouched phase is
    % perterbed.
    title({num2str(pv),num2str(20*log10(1./max(abs(yout))))})
    if(sum(ov) == N)
        disp('Total Loss (dB)');
        disp(20*log10(1./max(yout)))
        break;
    end
    pv(d1i) = pv(d1i) + .01;
    
    
    
    drawnow;
    

end




