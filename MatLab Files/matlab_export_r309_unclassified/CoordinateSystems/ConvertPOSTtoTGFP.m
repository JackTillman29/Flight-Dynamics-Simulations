clear variables;

% NOTE: all NASA POST output in feet, feet/sec & degrees. Do not convert
% to meters and radians until final output

% declare various conversions
deg2rad = pi / 180.0;
ft2m = 0.3048;

% declare constants consistent with NASA POST
equatorRadius = 2.0925741E7;
polarRadius = 2.085559E7;
omega = 7.29211E-5;  % rad/sec
omegaVector = [ 0 0 omega ]';

% compute Earth eccentricity
e = sqrt(equatorRadius.^2 - polarRadius.^2) / equatorRadius;

% get tab names from Excel spreadsheet
spreadsheet = 'AFcases2015Mar - part2.xlsx'
[status,sheets] = xlsfinfo(spreadsheet);

% loop on number of tabs
for sheetNum = 3:length(sheets)
    % load NASA POST data output
    [data,txt,raw]=xlsread(spreadsheet,sheetNum);
    
    % determine number of lines of data
    dataLines = length(data(:,1));
    
    % Use geodetic latitude, longitude & altitude of impact point to
    % define ECEF position of terminal point
    terminalLat = data(dataLines,18) * deg2rad;
    terminalLon = data(dataLines,19) * deg2rad;
    terminalAlt = data(dataLines,17);
    tempN = equatorRadius / sqrt(1.0 - (e*sin(terminalLat)).^2);
    terminalECEF = [ (tempN + terminalAlt)*cos(terminalLat)*cos(terminalLon) ...
        (tempN + terminalAlt)*cos(terminalLat)*sin(terminalLon) ...
        (tempN*(1-e.^2) + terminalAlt)*sin(terminalLat) ]';
    
    % Use geodetic latitude, longitude & 0 MSL of initial point to
    % define inertial Launcher coordinate system
    launcherLat = data(1,18) * deg2rad;
    %launcherLat = 0.0;
    launcherLon = data(1,19) * deg2rad;
    %launcherLon = 0.0;
    launcherAlt = 0.0;
    
    % define rotation matrix from ECI to launcher
    % NOTE: assumes azimuth of launcher c-sys equals zero
    TeciTOl = [ ...
        [ cos(launcherLat)*cos(launcherLon)   cos(launcherLat)*sin(launcherLon)  sin(launcherLat) ];
        [ -sin(launcherLon)                            cos(launcherLon)                     0     ];
        [ -sin(launcherLat)*cos(launcherLon) -sin(launcherLat)*sin(launcherLon)  cos(launcherLat) ] ];
    
    % define position of launcher c-sys relative to ECI
    positionL = [ (tempN + launcherAlt)*cos(launcherLat)*cos(launcherLon) ...
        (tempN + launcherAlt)*cos(launcherLat)*sin(launcherLon) ...
        (tempN*(1-e.^2) + launcherAlt)*sin(launcherLat) ]';
    
    outputData = [];
    for i = 1:dataLines
        % compute rotation matrix from ECI to ECEF
        time = data(i,1);
        TeciTOecef = [ ...
            [ cos(omega*time)  sin(omega*time) 0 ];
            [ -sin(omega*time) cos(omega*time) 0 ];
            [ 0                0               1 ] ];
        
        % compute rotation matrix from ECEF to E-N-U at terminal point
        geoLat = data(i,18) * deg2rad;
        geoLon = data(i,19) * deg2rad;
        TecefTOenu = [ ...
            [       -sin(geoLon)               cos(geoLon)              0      ];
            [ -cos(geoLon)*sin(geoLat)  -sin(geoLon)*sin(geoLat)   cos(geoLat) ];
            [  cos(geoLon)*cos(geoLat)   sin(geoLon)*cos(geoLat)   sin(geoLat) ] ];
        
        % get position in ECI frame
        positionECI = [ data(i,2) data(i,3) data(i,4) ]';
        
        % compute position of current point in ENU system
        positionENU = TecefTOenu.' \ (TeciTOecef.' \ positionECI - terminalECEF);
        
        % compute rotation matrix from ECI to ENU system
        TeciTOenu = TeciTOecef * TecefTOenu;
        
        % get velocity in ECI frame
        velocityECI = [ data(i,5) data(i,6) data(i,7) ]';
        
        % compute velocity of current point in ENU system
        velocityENU = TeciTOenu.' \ (velocityECI - cross((TeciTOecef.' * omegaVector), ...
            (TeciTOecef.' * terminalECEF + TeciTOenu.' \ positionENU)));
        
        % get Euler angles in launcher coordinate system
        roll = data(i,11) * deg2rad;
        yaw = data(i,12) * deg2rad;
        pitch = data(i,13) * deg2rad;
        
        % create body transformation matrix using 1-3-2 rotation (NASA POST)
        TlTOb = [ ...
            [ cos(yaw)*cos(pitch)  (cos(roll)*sin(yaw)*cos(pitch) + sin(roll)*sin(pitch)) ...
            (sin(roll)*sin(yaw)*cos(pitch) - cos(roll)*sin(pitch)) ];
            [ -sin(yaw)                cos(roll)*cos(yaw)      sin(roll)*cos(yaw) ];
            [ cos(yaw)*sin(pitch)  (cos(roll)*sin(yaw)*sin(pitch) - sin(roll)*cos(pitch)) ...
            (sin(roll)*sin(yaw)*sin(pitch) + cos(roll)*cos(pitch)) ] ];
        
        % create vectors in launcher c-sys for current orientation
        orientLx = TlTOb.' * [ 1 0 0 ]';
        orientLy = TlTOb.' * [ 0 1 0 ]';
        orientLz = TlTOb.' * [ 0 0 1 ]';
        
        % transform vectors into ENU system
        orientENUx = TecefTOenu.' \ ( TeciTOecef.' \ TeciTOl \ orientLx);
        orientENUy = TecefTOenu.' \ ( TeciTOecef.' \ TeciTOl \ orientLy);
        orientENUz = TecefTOenu.' \ ( TeciTOecef.' \ TeciTOl \ orientLz);
        
        % create rotation matrix in ENU system
        orientENU(1:3,1) = orientENUx;
        orientENU(1:3,2) = orientENUy;
        orientENU(1:3,3) = orientENUz;
        
        % NASA POST c-sys has body Z axis going out bottom;
        % to align with ESAMS, need to rotate about X-axis 180 degrees;
        orientENU = orientENU * [ [ 1 0 0 ]; [ 0 -1 0 ]; [ 0 0 -1 ] ];
        
        % extract Euler angles from rotation matrix
        yaw = atan2(orientENU(1,2),orientENU(1,1));
        pitch = asin(-orientENU(1,3));
        roll = atan2(orientENU(2,3),orientENU(3,3));
        
        % store new data
        outputData(i,1) = time;
        outputData(i,2) = positionENU(1) * ft2m;
        outputData(i,3) = positionENU(2) * ft2m;
        outputData(i,4) = positionENU(3) * ft2m;
        outputData(i,5) = velocityENU(1) * ft2m;
        outputData(i,6) = velocityENU(2) * ft2m;
        outputData(i,7) = velocityENU(3) * ft2m;
        outputData(i,8) = pitch;
        outputData(i,9) = roll;
        outputData(i,10) = yaw;
    end
    
    outputFile = strcat(sheets(sheetNum),'.txt');
    dlmwrite(char(outputFile),outputData,'delimiter','\t','precision',8);
end

