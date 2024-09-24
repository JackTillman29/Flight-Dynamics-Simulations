close all; clear all; clc

addpath('I:\MATLAB\Antenna\ARRAY_BLUE\') % latest phased array code
addpath('I:\MATLAB\Printing\')

c = 3e8;
r2d = 180/pi;
d2r = pi/180;

dbi2lin = @(x) 10.^(x/10);
lin2dbi = @(x) 10*log10(abs(x));
round2place = @(x,p) round(x.*10.^(p)) .* 10.^(-p);


arr = array('rectangular phased array');

rf = 10e9;
% rf = 600e6;
wavelen = c/rf;
spacefacx = 0.5;
spacefacy = 0.5;

efficiency = 1;

Nx = 10;
Ny = 10;

arr.CreateRectangularArray(wavelen,efficiency,spacefacx,spacefacy,Nx,Ny);
arr.recenterGeometry;
arr.compUniformWgt();

installYaw_deg   = 0;
installPitch_deg = 0;
installRoll_deg  = 0;
arr.changeInstNormal(installYaw_deg, installPitch_deg, installRoll_deg);

steerAz_deg = 0;
steerEl_deg = 0;
[azAnt,elAnt] = arr.xformAzEl(steerAz_deg * d2r, steerEl_deg * d2r);
arr.setSteer(azAnt,elAnt)



azext = 90*[-1 1];
elext = 90*[-1 1];
daz = 1;
del = 1;
% azext = 10*[-1 1];
% elext = 10*[-1 1];
% daz = 0.1;
% del = 0.1;
az = [azext(1):daz:azext(2)]; naz = length(az);
el = [elext(1):del:elext(2)]; nel = length(el);
[AZ,EL] = meshgrid(az,el);
% prod(size(AZ))
% return
g = arr.getGain(AZ*d2r,EL*d2r,arr.ENUM_MODE_AMPFILE);

g = abs(g).^2;
% figure; imagesc(az,el,10*log10(g))
figure; imagesc(g)
maxGain = max(max(10*log10(g)));
caxis(maxGain + [-50 0])


return
% write out antenna pattern for ESAMS
fid = fopen('tx_antenna.txt','w+');
fprintf(fid,'%s\n',num2str(naz));
fprintf(fid,'%s\n',num2str(nel));


for k = 1:nel
    fprintf(fid,'%s\n',num2str(el(k)));
end
for k = 1:naz
    
    fprintf(fid,'%s\n',num2str(az(k)));
    for j = 1:nel
        fprintf(fid,'%s\n',num2str(g(j,k)));
    end
    
end
fclose(fid)

