% simple script to draw spectral folding phenomena
% everything will be normalized to 1 = Fs, 1 = Max Amplitude

nTriangles = 2; % +/- 0Hz

% Analog spectra representation
analogFc = 190.1; % x Fs
analogBW = .2; % x Fs
analogLowAmp = 0.8;
analogHiAmp  = 1.0;

% generate "guide" sharktooth data
xData = (-nTriangles-0.5):0.5:(nTriangles+0.5);
yData = (mod(xData,1) == 0);

xAnalogLeft = [ ...
    -analogFc-analogBW/2 ...
    -analogFc-analogBW/2 ...
    -analogFc+analogBW/2 ...
    -analogFc+analogBW/2 ];
xAnalogRight = [ ...
     analogFc-analogBW/2 ...
     analogFc-analogBW/2 ...
     analogFc+analogBW/2 ...
     analogFc+analogBW/2];
 
 yAnalogLeft = [ ...
     0 analogHiAmp analogLowAmp 0 ];
 yAnalogRight = [ ...
     0 analogLowAmp analogHiAmp 0 ];
 
 % horizontal line


hf = figure(1);
hSharkTooth = plot(xData,yData);
%grid on;
set(gca,'NextPlot','add');
nReplications = round(analogFc);
for i = (-nReplications-nTriangles):(nReplications+nTriangles)
    if ( i == 0 )
        line(xAnalogRight-i,yAnalogRight,'Color','r', ...
            'LineStyle','-','LineWidth',2);
        line(xAnalogLeft-i, yAnalogLeft,'Color','g', ...
            'LineStyle','-','LineWidth',2);
    else
        line(xAnalogRight-i,yAnalogRight,'Color','r', ...
            'LineStyle',':','LineWidth',2);
        line(xAnalogLeft-i, yAnalogLeft,'Color','g', ...
            'LineStyle',':','LineWidth',2);
    end
end

set(gca,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1]);
set(gca,'XLim',[-nTriangles-.5 nTriangles+.5]);
grid on;
xlabel('Frequency (xFs)');
ylabel('Amplitude (Normalized)');

set(gca,'NextPlot','replace');