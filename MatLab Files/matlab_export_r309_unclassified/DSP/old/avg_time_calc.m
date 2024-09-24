close all; clear all; clc

fc = 1e6;

Fs = 100*fc; dt = 1/Fs;

T = 1e-3;
t = [0:dt:T-dt];
N = length(t);

a = linspace(1,0,N);
x = a.*cos(2*pi*fc.*t);
x = a;
% delta function
x = zeros(1,N);
x(randi(N,1)) = 1;
% pulse function
x = zeros(1,N);
stepWidth = 0.2*T
kWid = ceil(stepWidth / dt);
rndIdx = randi(N-kWid,1);
x(rndIdx:(rndIdx+kWid-1)) = 1;

xtd = abs(x).^2; % time density function
E = sum(xtd)*dt
xtd = xtd / E; % divide by total signal energy to normalize density function
% xtd = xtd ./ normFactor;

meanTint = cumsum(t.*xtd.*dt);     % <w>
ip_tsq   = meanTint.^2;            % <w>.^2
ip_t2    = cumsum(t.^2.*xtd.*dt);  % <w^2>  ip == inner product

% variance time == <w^2> - <w>^2
varTint = ip_t2 - ip_tsq;
% std time (RMS Time Duration)
stdTint = sqrt(varTint);

T
meanT = meanTint(end)
stdT = stdTint(end)


figure;
subplot(3,1,1)
plot(t,x); title(['signal, dur: ',num2str(T)])
subplot(3,1,2)
plot(t,meanTint); title(['mean time: ',num2str(meanT),' percent: ',num2str(meanT/T * 100)])

