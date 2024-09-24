a = array('a');


% (this,lambda,eff,spacefacx,spacefacy,Nx,Ny)

a.CreateRectangularArray(3e8/1e6, 1.0, 0.5, .5, 8, 1)
g1 = a.getAzCut((-90:1:90)*pi/180,0,a.ENUM_MODE_UNIFORM,'P');

a.CreateRectangularArray(3e8/100e6, 1.0, 0.5, .5, 8, 1)
g2 = a.getAzCut((-90:1:90)*pi/180,0,a.ENUM_MODE_UNIFORM,'P');

a.CreateRectangularArray(3e8/1e9, 1.0, 0.5, .5, 8, 1)
g3 = a.getAzCut((-90:1:90)*pi/180,0,a.ENUM_MODE_UNIFORM,'P');

figure;
hold on;
plot([-90:1:90],g1);
plot([-90:1:90],g2,'r');
plot([-90:1:90],g3,'bo');


ylim([-40 20]);
legend('100MHz','300MHz','900MHz');
