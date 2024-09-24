% close all; clear all; clc
clear all; clc

Fs = 40e3;
dt = 1/Fs;
% dt = 0.00001;

f = 100;
T = 0.05;
t = [0:dt:T-dt];

phs = 0.50;

% x1 = sawOsc(f,T,dt);
% x2 = sawOsc(f,T,dt,phs);

% dcy = 0.1;
% x1 = rectOsc(f,T,dcy,dt);
% x2 = rectOsc(f,T,dcy,dt,phs);

x1 = triOsc(f,T,dt);
x2 = triOsc(f,T,dt,phs);

% x3 = rectOsc(f,T,0.5,dt);
% x4 = rectOsc(f,T,0.5,dt,phs);

figure;
plot(t,x1)
hold on;
plot(t,x2,'r')
% xlim([0 1e-3])


return

%% sweep duty cycle
for dcy = linspace(0.1,0.9,300)
    x1 = rectOsc(f,T,dcy,dt);
    x2 = clockOsc(f,T,dcy,dt);
    if(~exist('hfig'))
        hfig = figure;
        hp1 = plot(t,x1);
        hold on;
        hp2 = plot(t,x2,'r');
    else
        set(hp1,'YData',x1)
        set(hp2,'YData',x2)
    end
    ylim([-2 2])
    drawnow;
end

%% sweep phase offset
dcy = 0.1;
for phs = linspace(0,1,3000)
%     x1 = rectOsc(f,T,dcy,dt);
%     x2 = rectOsc(f,T,dcy,dt,phs);
    
    x1 = triOsc(f,T,dt);
    x2 = triOsc(f,T,dt,phs);
    
    if(~exist('hfig'))
        hfig = figure;
        hp1 = plot(t,x1);
        hold on;
        hp2 = plot(t,x2,'r');
    else
        set(hp1,'YData',x1)
        set(hp2,'YData',x2)
    end
    ylim(1.1*[-1 1])
    drawnow;
end
