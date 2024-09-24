%Booz Allen Hamilton
%Jack Tillman

%This script stores the initial values and data for the 3DOF SimuLink model
%To modify parameters in the model, edit their numerical values here

%This script will use KMS (kilogram, meter, second) units
%All angle units can be input as degrees (a conversion to radians will
%occur at the end of the script!

%Clearing workspace and closing all possible outputs
clc;
clear all;
close all;

%Turning off warnings to suppress dialogue when parsing in data files
warning('off');

%%Simulation Time Values
initial_time = 0;
final_time = 100;
timestep = 0.001;

%%Atmospheric & Planetary Constants

%Gas Constant (Air)
rAir = 287.507;

%Cp/Cv for Air
gamma = 1.4;

%Spherical Radius of Earth
rEarth = 6371.0088 * 10^3;

%Elliptical Axes for Earth
aEarth = 6378.1370 * 10^3;
bEarth = 6356.7523 * 10^3;

%Earth rotation rate
omegaEarth = 7.2921159 * 10^-5;

%%Simulation Initial Conditions

%Initial mass & moments of inertia
m0 = 50;
Ixx = 10;
Iyy = 10;
Izz = 10;
Ixz = 1;

%Wing Area
wingArea = 50.0;

%Initial position in geodetic space
latitude0 = 45;
longitude0 = 45;
height0 = 0;

%Initial velocity in Wind Axes
u0Wind = 5;
v0Wind = 0;
w0Wind = 0;

%Initial position in Wind Axes
x0Wind = 0;
y0Wind = 0;
z0Wind = 0;

%Initial angular velocity in Wind Axes
p0Wind = 0.05;
q0Wind = 0;
r0Wind = 0.005;

%Initial angular position in Wind Axes
rollAngle0Wind = 0;
AoA0Wind = 0;
sideslipAngle0Wind = 0;

%%Force Coefficient & Altitude Values
    %These are read in through .csv files, so these need extra
    %consideration & computations. This correlates with more lines of code

%US Standard Atmosphere Data Tables & Interpolations
usAtmosphere = readtable("usAtmosphere.csv");

heightsList = usAtmosphere.Z;
densityList = usAtmosphere.x_;
temperatureList = usAtmosphere.T;

%Using Missile DatCom data for more accurate lift & drag coeff values
file_location = 'C:\Users\641893\BAH Files\Missile DATCOM\DATCOM Testing\Trimmed Coeff Testing\';

%Specifying the file name (This doesn't change & is set by MDC!)
file_name = 'for042.csv';

%Joining the two strings to specify the complete file path
file_path = strcat(file_location, file_name);

%Loading in the data as a table:
data_table = rmmissing(readtable(file_path, 'VariableNamingRule','preserve'));

%storing missile datcom data as arrays - these arrays are used to
%interpolate in the SimuLink model
alpha_MDC = table2array(data_table(:, '''ALPHA'''));
mach_MDC = table2array(data_table(:, '''MACH'''));
altitude_MDC = table2array(data_table(:, '''ALT'''));
lift_coeffs_MDC = table2array(data_table(:, '''TRIM_CL'''));
drag_coeffs_MDC = table2array(data_table(:, "'TRIM_CD'"));

lift_coeff_table = table(altitude_MDC, mach_MDC, alpha_MDC, lift_coeffs_MDC);
drag_coeff_table = table(altitude_MDC, mach_MDC, alpha_MDC, drag_coeffs_MDC);

lift_coeffs_array = [altitude_MDC,mach_MDC,alpha_MDC,lift_coeffs_MDC];
drag_coeffs_array = [altitude_MDC,mach_MDC,alpha_MDC,drag_coeffs_MDC];

%Now to format this data such that it's usable for SimuLink

numAltitudes = 21;
numAoAs = 13;
numMachs = 20;

%initializing our looping variables
i = 1; %i loops through each row
j = 1; %j loops through each column
k = 1; %k loops through each stack of rectangles

%pre - initializing our 3D array
cL_3DArray = zeros(numAoAs, numMachs, numAltitudes);
cD_3DArray = zeros(numAoAs, numMachs, numAltitudes);

%initializing looping values
loopAlt = 0;
loopAngle = -30;
loopMach = 0.5;
loopVariable = 1;

while loopAlt <= 60000
    loopMach = 0.5;
    while loopMach <= 10
        loopAngle = -30;
        while loopAngle <= 30.0
            cL_3DArray(i, j, k) = lift_coeffs_array(loopVariable, 4);
            loopVariable = loopVariable + 1;
            i = i + 1;
            loopAngle = loopAngle + 5.0;
        end
        i = 1;
        j = j + 1;
        loopMach = loopMach + 0.5;
    end
    j = 1;
    k = k + 1;
    if loopAlt < 10000
        loopAlt = loopAlt + 1000;
    elseif loopAlt >= 10000
        loopAlt = loopAlt + 5000;
    end
end

%initializing looping values
loopAlt = 0;
loopAngle = -30;
loopMach = 0.5;
loopVariable = 1;
i = 1;
j = 1;
k = 1;


while loopAlt <= 60000
    loopMach = 0.5;
    while loopMach <= 10
        loopAngle = -30;
        while loopAngle <= 30.0
            cD_3DArray(i, j, k) = drag_coeffs_array(loopVariable, 4);
            loopVariable = loopVariable + 1;
            i = i + 1;
            loopAngle = loopAngle + 5.0;
        end
        i = 1;
        j = j + 1;
        loopMach = loopMach + 0.5;
    end
    j = 1;
    k = k + 1;
    if loopAlt < 10000
        loopAlt = loopAlt + 1000;
    elseif loopAlt >= 10000
        loopAlt = loopAlt + 5000;
    end
end

machArray = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, ...
    8.5, 9.0, 9.5, 10.0];
altitudeArray = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, ...
    25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000];
alphaArray = [-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30];

%Performing a similar procedure to store the optimal angles of attack

%Initializing arrays that'll be used to store data
LD_Ratio = zeros(numAoAs, numMachs, numAltitudes);
optimalAoA = zeros(numMachs, numAltitudes);

%loop variable resetting
loopAlt = 0;
loopAngle = -30;
loopMach = 0.5;
loopVariable = 1;
i = 1;
j = 1;
k = 1;

while loopAlt <= 60000
    loopMach = 0.5;
    while loopMach <= 10
        loopAngle = -30;
        while loopAngle <= 30.0
            LD_Ratio(i, j, k) = lift_coeffs_array(loopVariable, 4) / drag_coeffs_array(loopVariable, 4);
            loopVariable = loopVariable + 1;
            i = i + 1;
            loopAngle = loopAngle + 5.0;
        end
        i = 1;
        j = j + 1;
        loopMach = loopMach + 0.5;
    end
    j = 1;
    k = k + 1;
    if loopAlt < 10000
        loopAlt = loopAlt + 1000;
    elseif loopAlt >= 10000
        loopAlt = loopAlt + 5000;
    end
end

%Now we just need to maximize each column in the LD_Ratio array, locate the corresponding
%index, compute the corresponding AoA, & then store these values in the
%optimalAoA array

%loop variable resetting
loopAlt = 0;
loopAngle = -30;
loopMach = 0.5;
loopVariable = 1;
i = 1;
j = 1;
k = 1;

while loopAlt <= 60000
    loopMach = 0.5;
    while loopMach <= 10
        loopAngle = -30;
        LDRatioTemp = LD_Ratio(:, j, k);
        maxLDRatioTemp = max(LDRatioTemp);
        maxLDRatio_Index = find(LDRatioTemp == maxLDRatioTemp);
        maxAoA = -30 + 5 * (maxLDRatio_Index - 1);
        optimalAoA(j, k) = maxAoA;
        j = j + 1;
        loopMach = loopMach + 0.5;
    end
    j = 1;
    k = k + 1;
    if loopAlt < 10000
        loopAlt = loopAlt + 1000;
    elseif loopAlt >= 10000
        loopAlt = loopAlt + 5000;
    end
end


%%Unit Conversions

%Degrees to Radians Conversion
function radians = deg2rad(x)
    radians = x * pi/180;
end

%Initial latitude & longitude conversions
latitude0 = deg2rad(latitude0);
longitude0 = deg2rad(longitude0);

%Initial angular position conversions (wind axes)
rollAngle0Wind = deg2rad(rollAngle0Wind);
AoA0Wind = deg2rad(AoA0Wind);
sideslipAngle0Wind = deg2rad(sideslipAngle0Wind);

%Initial angular velocity conversions (wind axes)
p0Wind = deg2rad(p0Wind);
q0Wind = deg2rad(q0Wind);
r0Wind = deg2rad(r0Wind);

%Printing message to screen if successful
successMessage = "All constants & values loaded in successfully!" + newline;
disp(successMessage);