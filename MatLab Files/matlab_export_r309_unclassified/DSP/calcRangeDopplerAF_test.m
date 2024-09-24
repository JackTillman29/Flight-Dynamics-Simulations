close all; clear all; clc

pf = @(x) 10*log10(abs(x));
c = 3e8;

Fs = 10e6; dt = 1/Fs;
IF = 0;
bw = 1.2e6;
dur = 512e-6;

t = [0:dt:dur-dt];
N = length(t);

f = linspace(-bw/2,bw/2,N);
x = exp(1i*2*pi*cumsum(f).*dt);

% figure; plot(real(x))

dopFac = 0.25;
dopExtent = (dopFac*Fs)*[-1 1];
nDop = 1000;
tic
% [af,rd,dopArray] = calcRangeDopplerAF(x,x,dt);
[af,rd,dopArray] = calcRangeDopplerAF(x,x,dt,dopExtent,nDop);
toc
% return

clims = [-10 0];

figure; add_print_callbacks;
imagesc(rd/c, dopArray, pf(abs(af)).')
imagesc(pf(abs(af)).')
xlabel('t_d [s]')
ylabel('Doppler [Hz]')
caxis(clims)
drawnow;

% return

tic
nDop = length(dopArray);
kdop = 0;
for dop = dopArray
    kdop = kdop + 1;
    disp(['kdop = ',num2str(kdop),' / ',num2str(nDop)])
    dopsig = exp(-1i*2*pi*dop.*t);
    tmp = PulseCompressionXCORRsStruct(x.*dopsig, x);
    if(~exist('af2'))
%         af2 = zeros(length(tmp.signalCrossCorTrim),nDop);
        af2 = zeros(length(tmp.signalCrossCor),nDop);
    end
    
%     af2(:,kdop) = tmp.signalCrossCorTrim.';
    af2(:,kdop) = tmp.signalCrossCor.';
    
end
toc


figure; add_print_callbacks;
subplot(2,1,1)
plot(rd/c*1e6,pf(af(:,500)));
hold on;
plot(rd/c*1e6,pf(af2(:,500)));
subplot(2,1,2)
plot(pf(af(500,:)));
hold on;
% plot(pf(af2(500,:)));
xlabel('t_d [\mus]')

figure; add_print_callbacks;
imagesc(1:size(af2,1), dopArray, pf(abs(af2)).')
% imagesc(pf(abs(af2)).')
xlabel('t_d [s]')
ylabel('Doppler [Hz]')
caxis(clims)



