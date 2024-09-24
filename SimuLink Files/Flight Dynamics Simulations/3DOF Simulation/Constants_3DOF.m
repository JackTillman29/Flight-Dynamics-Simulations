%Booz Allen Hamilton
%Jack Tillman

%This script stores the initial values and data for the 3DOF SimuLink model
%To modify parameters in the model, edit their numerical values here

%This script will use KMS (kilogram, meter, second) units

%Clearing workspace and closing all possible outputs
clc;
clear all;
close all;

%Simulation Time Values
initial_time = 0;
final_time = 700;
timestep = 0.001;

%Universal Constants & Values

rAir = 287.507; %Gas Constant of Air
gamma = 1.4; %Cp/Cv for air

mu_earth = 3.986004418 * 10^14;

r_earth = 6371.0088 * 10^3;

omega_earth = 7.2921159 * 10^-5; %Earth rotation rate

%Sea Level Atmospheric Conditions
rho_sl = 1.225;
temperature_sl = 288.15;
pressure_sl = 10.13 * 10^4;

%Initial Conditions

%Input Angles in degrees, these may be converted into radians if necessary.
%Input any values in terms of their KMS (kilogram, meter, seconds) units. 

wing_area = 0.079813;
mass = 161.5;

speed_0 = 1500;

height_0 = 35 * 10^3;
range_to_go_0 = 11112 * 10^3;
range_0 = 0;

latitude_0 = 38.99;
longitude_0 = -72.42;

heading_angle_0 = 20;
flight_path_angle_0 = 3;


%Degrees to radians conversions
latitude_0 = latitude_0 * (pi/180);
longitude_0 = longitude_0 * (pi/180);

heading_angle_0 = heading_angle_0 * (pi/180);
flight_path_angle_0 = flight_path_angle_0 * (pi / 180);

%US Standard Atmosphere Data Tables & Interpolations
usAtmosphere = readtable("usAtmosphere.csv");

heightsList = usAtmosphere.Z;
rhoList = usAtmosphere.x_;
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

%Formatting data such that it's valid for SimuLink lookup table

%Now let's format this data into a 3D array, this will allow SimuLink to
%interpolate between points in the array

%It will have the mach number on the x-axis, alpha on the y and lastly
%altitude on the z

%as such, we can write our array in this manner. since we know the total
%number of each different variable, we can easily use for loops

numAltitudes = 21;
numAoAs = 13;
numMachs = 20;

%Our matrix will look like stacks of rectangles, 21 to be exact, with each
%rectangle acting as a 13x20 rectangle

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

%We can analyze the lift and drag coefficients to tabulate optimal angles
%of attack for different flight conditions. Doing so, we can interpolate
%between values to ensure that the FVD simulation flies the aircraft at the
%most optimal AoA

%To do this, we need to create a table of AoA values for each flight
%parameter. To find the optimal AoA, we need to consider each situation of
%steady level flight -- an initial mach and altitude

%We just need to optimize each column in our 3D Arrays! This reduces our
%number of rows down from 13 (13 AoAs) to just 1 (AoA_Opt). We can use this
%information to initialize our new 3D matrices:

LD_Ratio = zeros(numAoAs, numMachs, numAltitudes);
optimalAoA = zeros(numMachs, numAltitudes);

%We can loop through this using a similar algorithm as before:

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

%Now we just need to maximize each column in the LD_Ratio column & find
%what AoA that corresponds to -- sort through each column, maximize, locate
%the index, compute the AoA that index corresponds with, store in
%optimalAoA array

%looping through every column for each altitude: 

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

%Now that we tabulated our maximum AoAs, we can use this to dynamically
%optimize the angle of attack in our simulink model!

