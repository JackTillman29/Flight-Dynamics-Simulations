%Jack Tillman
%Booz Allen Hamilton

%This script parses in data, analyzes and & plots multiple different
%figures. Since this script is parsing in data from a 6DOF flight
%simulation, there is a large amount of data files & plots to compute. 

%Things to note when running this script:

    %All .mat data files must be saved as arrays in the SimuLink model --
    %timeseries data Will Throw Errors!

    %Run this file immediately after the conclusion of a SimuLink
    %simulation, this will move the data files & plots to a different
    %directory, essentially cleaning this one in the process

    %This script will only work on simulation files that were created
    %successfully - no errors must have occurred in the simulation for this
    %script to function properly

    %Triple check time values to make sure they match the ones used in the
    %SimuLink simulation, if they dont match exactly -> errors!

    %All angular values are in degrees - this makes it a bit easier to read
    %and interpret!

%Clearing workspace and closing all possible outputs
clc;
clear all;
close all;

%Simulation Time Values
initial_time = 0;
final_time = 150;
timestep = 0.001;

%Data Parsing mode - 0 means files were already parsed & moved, 1 means the
%files must be moved... use 0 if this script has already been ran and 1 if
%running for the first time following a SimuLink simulation
dataParsingMode = 1;

%To parse in our data files, we'll specify their file path and take
%advantage of the fact that all files start with the same string: 6DOF_...

%We can then specify the names of each variable we're considering and
%concatenate them with the file path string to load in each file!

%Initial Path location
init_data_folder = "C:\Users\641893\BAH Files\MatLab & SimuLink Files\SimuLink Files\Flight Dynamics Simulations\6DOF Simulation\6DOF_";

%Final path location
final_data_folder = "C:\Users\641893\BAH Files\MatLab & SimuLink Files\SimuLink Files\Flight Dynamics Simulations\6DOF Simulation\6DOF_Data\";

%Array of the variable names (unfortunately it's a ton..)
data_names_array = ["AoA", "bankAngle","density", "drag", "dragCoeff", "earthRotation", "ecefX", "ecefY", "ecefZ",...
    "flightpathAngle", "geodeticRadius", "gravity", "inertialPitch", "inertialRoll", ...
    "inertialX", "inertialY", "inertialYaw", "inertialZ", "latitude", "lift", "liftCoeff", "longitude",...
    "mach", "mass", "nedX", "nedY", "nedZ", "pBody", "phiBody", "pitchMoment", "pitchSurfaceDeflection", "psiBody",...
    "qBody", "rBody", "rollMoment", "rollSurfaceDeflection", "sideslip", "thetaBody", ...
    "thrust", "thrustAoA", "thrustSideslip", "time", "uBody", "vBody", "velocity", "wBody",...
    "xForce", "yawMoment", "yawSurfaceDeflection", "yForce", "zForce"];

%Data file extension
data_extension = ".mat";

%Let's now make the directory where we'll store the .pngs of the plots we
%create

plotDirectory = final_data_folder + "Plots & Figures\";

%We can now initialize the array that'll be used to store all of our data

%Computing size of array
numTimeSteps = 1 + ((final_time - initial_time) / timestep);
numDataNames = length(data_names_array);

%Initializing array storing all of our data
data_array = zeros(numDataNames, numTimeSteps);

%We can now store all of our data into this array (via a loop)

%Looping through our data files

%Mode 0

%Loop indices
i = 1;
final_index = length(data_names_array);

if dataParsingMode == 0
    for i = 1:final_index

        %Loading in our data files
        file_name = strcat(final_data_folder + "6DOF_", data_names_array(i), data_extension);
        temp_data10 = struct2cell( load(file_name) );
        temp_data20 = temp_data10{1,1};

        %Writing data files to larger array
        data_array(i, :) = temp_data20(2, :);
    end
end

%Mode 1

%Loop indices
i = 1;
final_index = length(data_names_array);

if dataParsingMode == 1
    for i = 1:final_index

        %Loading in our data files
        file_name = strcat(init_data_folder, data_names_array(i), data_extension);
        temp_data11 = struct2cell( load(file_name) );
        temp_data21 = temp_data11{1,1};

        %Writing data files to larger array
        data_array(i, :) = temp_data21(2, :);

        %Moving data file to new location
        movefile(file_name, final_data_folder);

    end
    mkdir(plotDirectory);
end

%Printing message to screen if data parsing was successful
successMessage = "Data was parsed successfully!" + newline;
disp(successMessage);

%Now that we have all of our data stored, we can plot our values as we'd
%like

%Initial indexing
radiusIndex = find(data_names_array == "geodeticRadius");
timeIndex = find(data_names_array == "time");
inertialZIndex = find(data_names_array == "inertialZ");

%Let's first find the time at which the height of our vehicle is equal to
%0, this allows us to place labels on our plots at this point in time
time0 = interp1(data_array(inertialZIndex, :), data_array(timeIndex, :), 0);

%Since we need to create a bunch of plots, we can write functions to make
%our lives much easier!

function plotter = plottingFunction(data_array, xIndex, yIndex, xLabel, yLabel, plotTitle, plotLegend, fileLocation, height0Label, time0)
    figure

    plot(data_array(xIndex, :), data_array(yIndex, :), 'b');

    if height0Label == 1
        xline(time0, '--k', 'label', 'Height = 0 meters', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','middle');
    end

    xlabel(xLabel);
    ylabel(yLabel);
    grid on;
    title(plotTitle);
    legend(plotLegend, 'Location', 'best');
    plotFileName = fileLocation + "6DOF " + plotTitle + ".png";
    saveas(gcf, plotFileName);
end

%%Plotting the z axis height in the wind frame with respect to time
plottingFunction(data_array, timeIndex, inertialZIndex, "Time (s)", "Height (m)", "Height versus Time (W.A.)", "Height", plotDirectory, 1, time0);

%%Plotting the velocity with respect to time

%First grabbing the velocity index
velocityIndex = find(data_names_array == "velocity");

%Plotting
plottingFunction(data_array, timeIndex, velocityIndex, "Time (s)", "Velocity (m/s)", "Velocity versus Time", "Velocity", plotDirectory, 1, time0);

%%Let's now plot the geodetic radius with respect to time
plottingFunction(data_array, timeIndex, radiusIndex, "Time (s)", "Geodetic Radius (m)", "Geo. Radius versus Time", "Geodetic Radius", plotDirectory, 1, time0);

%%Let's plot our angle of attack as a function of mach number

%Grabbing the correct indices
aoaIndex = find(data_names_array =="AoA");
machIndex = find(data_names_array == "mach");

%Plotting
plottingFunction(data_array, machIndex, aoaIndex, "Mach Number", "Angle of Attack (Degrees)", "AoA versus Mach", "Angle of Attack", plotDirectory, 0, time0);

%%Let's now focus on plotting our range versus height, let's look at this
%%in the wind axis for now

%Indices
inertialXIndex = find(data_names_array == "inertialX");

%Plotting
plottingFunction(data_array, inertialXIndex, inertialZIndex, "Range (m)", "Height (m)", "Height versus Range (W.A.)", "Height", plotDirectory, 0, time0);

%%We can now plot the variation in our flight path angle with time

%Index!
flightPathAngleIndex = find(data_names_array == "flightpathAngle");

%Plotting
plottingFunction(data_array, timeIndex, flightPathAngleIndex, "Time (s)", "Flight Path Angle (Degrees)", "Flight Path Angle versus Time", "Flight Path Angle", plotDirectory, 0, time0);

%%We'll now also plot our bank angle with respect to time

%Index
bankAngleIndex = find(data_names_array == "bankAngle");

%Plotting
plottingFunction(data_array, timeIndex, bankAngleIndex, "Time (s)", "Bank Angle (Degrees)", "Bank Angle versus Time", "Bank Angle", plotDirectory, 0, time0);

%Plotting Earth's rotation rate with time - this should be constant with
%respect to time!

%Indices!
earthRotationIndex = find(data_names_array == "earthRotation");

%Plotting
plottingFunction(data_array, timeIndex, earthRotationIndex, "Time (s)", "Earth Rotation Rate (deg/s)", "Earth's Rotation versus Time", "Earth's Rotation", plotDirectory, 0, time0);

%Let's now plot our LLA Coordinates to see how this flight vehicle moved
%through Earth's atmosphere

%First, indices!
latitude_index = find(data_names_array == "latitude");
longitude_index = find(data_names_array == "longitude");
geoRadius_index = find(data_names_array == "geodeticRadius");

%plotting
uif = uifigure;
g = geoglobe(uif);
geoplot3(g, data_array(latitude_index, :), data_array(longitude_index, :), data_array(geoRadius_index, :), 'color', 'red', 'linewidth', 1.5);
