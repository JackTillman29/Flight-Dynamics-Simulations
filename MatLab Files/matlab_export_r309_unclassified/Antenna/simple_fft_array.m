close all;
clear;
clc;

nR = 21;
nC = 21;

g  = ones(nR,nC);
[xk,yk] = meshgrid(-(nR-1)/2:(nR-1)/2,-(nC-1)/2:(nC-1)/2);


% compute phasing
az = 0*pi/180;
el = 0*pi/180;

% phase noise
nz.phase_mag = 45*pi/180;
nz.amp_mag = 0.6; % percent
nz.phase = exp(nz.phase_mag*1i*2*(rand(nR,nC)-0.5));
nz.amp = 1.0 - nz.amp_mag*(rand(nR,nC));
rk = sqrt(xk.^2 + yk.^2);
rk = rk ./ max(rk(:));
rk = 1.0 - rk;
wgt = sin(pi/2*rk).^2;




% cos taper
%a = (cos(rg).^(1))';
% uniform taper
a = nz.amp .* wgt;

% broader spacing
if(0)
    a(2:2:end,:)=0;
    a(:,2:2:end)=0;
    uscale = 2;
    vscale = 2;
else
    uscale = 1;
    vscale = 1;
end

% compute peak amplitude (via constructive interference)
apeak = sum(a(:));

NFFT = 2048;
windowSize = NFFT/uscale;

pa = exp(1i*pi*sin(az) * yk/vscale);
pe = exp(1i*pi*sin(el) * xk/uscale);
ps = (pa.*pe).';

fg = fftshift(fft2(a.*g.*ps.*nz.phase,NFFT,NFFT));

u = linspace(-1,1,NFFT);
v = linspace(-1,1,NFFT);

% cut data
ucut = u(1:uscale:end);
vcut = v(1:vscale:end);
% how many copies of pattern present?
nCopiesX = uscale-1;
nCopiesY = vscale-1;
datSizeX = NFFT / uscale;
datSizeY = NFFT / vscale;
if(nCopiesX > 0)
    startu = round(NFFT / (nCopiesX * 4) + 1);
    startv = round(NFFT / (nCopiesY * 4) + 1);
    stopu = round(startu - 1 + datSizeX);
    stopv = round(startv - 1 + datSizeY);
else
    startu = 1;
    startv = 1;
    stopu = NFFT;
    stopv = NFFT;
end

fgcut = fg(startv:stopv,startu:stopu)./apeak;

maxfg = max(abs(fgcut(:)));

% compute w
[vcutg,ucutg]=meshgrid(vcut,ucut);
wcutg = sqrt(1.0 - ucutg.^2 - vcutg.^2);
wrcutg = real(wcutg);

% compute element pattern
angofb = acos(wrcutg);
elem_pat = cos(angofb);
fgcut_e = fgcut .* elem_pat;

figure;
subplot(1,2,1);
hi=imagesc(vcut,ucut,20*log10(abs(fgcut_e)));
set(gca,'YDir','normal');
caxis([-70 0]);
xlabel('U');
ylabel('V');
title('Sin Space');
axis equal;
axis tight;
colorbar;



% Warp to az el space
[el,az]=meshgrid(-90:.1:90,-90:.1:90);
el = el * pi/180;
az = az * pi/180;
[wlook,vlook,ulook]=sph2cart(az,el,1.0);
%[vin,uin]=meshgrid(linspace(-1,1,1024),linspace(-1,1,1024));
% 
gazel = interp2(vcutg,ucutg,abs(fgcut_e),vlook,ulook)';
% 

subplot(1,2,2);
hi=imagesc(180/pi*el(1,:),180/pi*az(:,1),20*log10(abs(gazel)));
set(gca,'YDir','normal');
caxis([-70 0]);
xlabel('Az (deg)');
ylabel('El (deg)');
title('Az/El Space');
axis equal;
axis tight;
colorbar;

try
    add_print_callbacks;
catch
    try
        addpath('z:\sawmillkd\MATLAB\Printing');
        add_print_callbacks;
    catch
        disp('could not load print callbacks!');
    end
end