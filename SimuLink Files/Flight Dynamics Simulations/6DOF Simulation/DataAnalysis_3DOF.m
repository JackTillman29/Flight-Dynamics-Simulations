%Jack Tillman
%Booz Allen Hamilton 

clc
clear all;
close all;

%This script parses in data collected from a SimuLink simulation and plots
%as needed. Data analysis will also be included as well.

%Parsing in the data files

%Path location
data_folder = "C:\Users\641893\BAH Files\MatLab & SimuLink Files\SimuLink Files\Flight Dynamics Simulations\3DOF Simulation\3DOF_Data\3DOF_";

%Data names
data_names_array = ["angle_of_attack" "bank_angle" "density" "drag" "drag_coeff" "flight_path_angle" "gravity" "heading_angle" "height" "latitude" "lift"...
    "lift_coeff" "longitude" "mach" "radius" "range_to_go" "speed" "range"];

%Data file extension
data_extension = ".mat";

%We can now loop through all of our data files and store them in a single
%large array. Doing so, we'll save all of the time and data values for each
%variabe we're considering. As such, each variable will hold 2 complete
%rows in our larger array

%Simulation Time Values
initial_time = 0;
final_time = 700;
timestep = 0.001;

numTimeSteps = 1 + ((final_time - initial_time) / timestep);
numDataNames = 1 + length(data_names_array);

%Data array containing all data (note, 18 indices is for the 17 variables +
%1 for the time
data_array = zeros(numDataNames, numTimeSteps);

%Writing the time to the first row in the array
data_array(1,:) = linspace(0, final_time, numTimeSteps);

%Loop indices
i = 1;
final_index = length(data_names_array);

%Looping through our data files
for i = 1:final_index

    %Loading in our data files
    file_name = strcat(data_folder, data_names_array(i), data_extension);
    temp_data1 = struct2cell( load(file_name) );
    temp_data2 = temp_data1{1,1};

    %Writing data files to larger array
    data_array(i + 1, :) = temp_data2(2, :);

end

%Now that we've stored everything to our list, we can plot our data.
%Finally. 

%We'll start with plotting the speed, to do that, we'll need to find what
%index is associated with that variable. We'll be careful to add 1 to
%account for the index associated with the time

speed_index = find(data_names_array == "speed") + 1;

%The index associated for the time is easy, it's just 1
time_index = 1;

%Finding height index:
height_index = find(data_names_array == "height") + 1;

%First, let's find the time that the height is equal to 0:
time0 = interp1(data_array(height_index, :), data_array(time_index, :), 0);

%We can now actually plot! We're now plotting a speed versus time plot
figure;
plot(data_array(time_index, :), data_array(speed_index, :), 'b');
xlabel('Time (s)');
ylabel('Speed (m/s)');
xlim([initial_time, final_time]);
xline(time0, '--k', 'label', 'Height = 0 meters', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','middle');
grid on;
title('Speed of Hypersonic Re-entry Vehicle (3DOF)')
legend(data_names_array(speed_index -1), 'Location','southwest');

%Let's now plot a speed versus height to see how they compare

%Plotting
figure;
plot(data_array(height_index, :), data_array(speed_index, :), 'b');
xlabel('Height (m)');
ylabel('Speed (m/s)');
xlim([-2000, max(data_array(height_index, :)) + 2000]);
set(gca,'xdir','reverse');
grid on;
title('Speed of Hypersonic Re-entry Vehicle (3DOF)')
legend(data_names_array(speed_index -1), 'Location','southwest');

%Plotting a height versus time graph
figure;
plot(data_array(time_index, :), data_array(height_index, :), 'b');
xlabel('Time (s)');
ylabel('Height (m)');
ylim([-2000, max(data_array(height_index, :)) + 2000]);
xlim([initial_time, final_time]);
grid on;
title('Height for a HRV');
legend(data_names_array(height_index - 1), 'Location', 'southwest');

%Let's now plot the drag and lift coeffs as functions of time, we'll do the
%same for height as well

drag_coeff_index = 1 + find(data_names_array == "drag_coeff");
lift_coeff_index = 1 + find(data_names_array == "lift_coeff");

figure;
plot(data_array(time_index, :), data_array(lift_coeff_index, :), 'b');
hold on;
plot(data_array(time_index, :), data_array(drag_coeff_index, :), 'Color', '#EDB120');
xline(time0, '--k', 'label', 'Height = 0 meters', 'LabelHorizontalAlignment','left');
xlabel("Time (s)");
ylabel("Aero Coefficients")
xlim([initial_time, final_time]);
grid on;
title('C_{l} & C_{d} of Hypersonic Re-entry Vehicle (3DOF)');
legend({data_names_array(lift_coeff_index - 1), data_names_array(drag_coeff_index - 1)}, 'Interpreter', 'none', 'Location', 'northwest');
hold off;

%Plotting these versus height

figure;
plot(data_array(height_index, :), data_array(lift_coeff_index, :), 'b');
hold on;
plot(data_array(height_index, :), data_array(drag_coeff_index, :), 'Color', '#EDB120');
xlabel("Height (m)");
ylabel("Aero Coefficients");
xlim([-2000, max(data_array(height_index, :)) + 2000]);
set(gca,'xdir','reverse');
grid on;
title('C_{l} & C_{d} of Hypersonic Re-entry Vehicle (3DOF)');
legend({data_names_array(lift_coeff_index - 1), data_names_array(drag_coeff_index - 1)}, 'Interpreter', 'none', 'Location', 'northwest');
hold off;

%Now plotting lift and drag versus time and then height

drag_index = 1 + find(data_names_array == "drag");
lift_index = 1 + find(data_names_array == "lift");

%time
figure;
plot(data_array(time_index, :), data_array(lift_index, :), 'b');
hold on;
plot(data_array(time_index, :), data_array(drag_index, :), 'Color', '#EDB120');
xline(time0, '--k', 'label', 'Height = 0 meters', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
xlabel("Time (s)");
ylabel("Aero Forces (N)")
xlim([initial_time, final_time]);
grid on;
title('Lift & Drag of Hypersonic Re-entry Vehicle (3DOF)');
legend({data_names_array(lift_index - 1), data_names_array(drag_index - 1)}, 'Interpreter', 'none', 'Location', 'northwest');
hold off;

%height
figure;
plot(data_array(height_index, :), data_array(lift_index, :), 'b');
hold on;
plot(data_array(height_index, :), data_array(drag_index, :), 'Color', '#EDB120');
xlabel("Height (m)");
ylabel("Aero Forces (N)");
xlim([-2000, max(data_array(height_index, :)) + 2000]);
set(gca, 'xdir', 'reverse');
grid on;
title('Lift & Drag of Hypersonic Re-entry Vehicle (3DOF)');
legend({data_names_array(lift_index - 1), data_names_array(drag_index - 1)}, 'Interpreter', 'none', 'Location', 'northwest');
hold off;

%now let's plot our angular values 

%starting with our latitude and longitude:

latitude_index = 1 + find(data_names_array == "latitude");
longitude_index = 1 + find(data_names_array == "longitude");

%plotting
uif = uifigure;
g = geoglobe(uif);
geoplot3(g, data_array(latitude_index, :), data_array(longitude_index, :), data_array(height_index, :), 'color', 'red', 'linewidth', 1.5);

%now, we can plot our flight path angle and heading angle, these can be
%done as a function of height as well as time

%since these are angles, we can plot these on polar plots to see if they
%are cyclic

%indices
flight_path_index = 1 + find(data_names_array == "flight_path_angle");
heading_index = 1 + find(data_names_array == "heading_angle");

%plotting with respect to time first
figure;
plot(data_array(time_index, :), data_array(flight_path_index, :), 'color', 'b')
hold on;
plot(data_array(time_index, :), data_array(heading_index, :), 'color', '#EDB120');
title('Flight Path & Heading Angles wrt to Time (3DOF)');
xlabel('Time (s)')
ylabel('Angles (degrees)');
xlim([initial_time, final_time]);
xline(time0, '--k', 'label', 'Height = 0 meters', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
grid on;
legend({data_names_array(flight_path_index - 1), data_names_array(heading_index - 1)}, 'Interpreter','none', 'Location','northwest');
hold off;

%plotting with respect to height
figure;
plot(data_array(height_index, :), data_array(flight_path_index, :), 'color', 'b')
hold on;
plot(data_array(height_index, :), data_array(heading_index, :), 'color', '#EDB120');
title('Flight Path & Heading Angles wrt to Height (3DOF)');
grid on;
xlim([-2000, max(data_array(height_index, :)) + 2000]);
set(gca, 'xdir', 'reverse');
xlabel('Height (m)')
ylabel('Angles (degrees)');
legend({data_names_array(flight_path_index - 1), data_names_array(heading_index - 1)}, 'Interpreter','none', 'Location','northwest');
hold off;

%plotting height versus range to go

%first, indices
range_to_go_index = 1 + find(data_names_array == "range_to_go");

%plotting
figure;
plot(data_array(range_to_go_index, :), data_array(height_index, :), 'b');
xlabel('Range To Go (m)');
ylabel('Height (m)');
grid on;
ylim([-2000, max(data_array(height_index, :)) + 2000]);
set(gca, 'xdir', 'reverse');
title('Height v Range To Go of a HRV (3DOF)');
legend(data_names_array(height_index -1), 'Location','southwest');

%finding range index
range_index = 1 + find(data_names_array == "range");
figure;
plot(data_array(range_index, :), data_array(height_index, :), 'b');
xlabel('Range (m)');
ylabel('Height (m)');
grid on;
ylim([-2000, max(data_array(height_index, :)) + 2000]);
title('Height v Range of a HRV (3DOF)');
legend(data_names_array(height_index -1), 'Location','southwest');

%Finding optimal alpha to maximize range -- max of L/D 

LD_Ratio = data_array(lift_index, :) / data_array(drag_index, :);
alphaOpt = (180 / pi) * atan(1/max(LD_Ratio));

%Plotting angle of attack versus height and time
alpha_index = 1 + find(data_names_array == "angle_of_attack");

figure;
plot(data_array(time_index, :), data_array(alpha_index, :), 'b');
xlabel('Time (s)');
ylabel('Speed (m/s)');
xlim([initial_time, final_time]);
ylim([5, max(data_array(alpha_index, :)) + 5]);
xline(time0, '--k', 'label', 'Height = 0 meters', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top');
grid on;
title('AoA of a Hypersonic Re-entry Vehicle (3DOF)')
legend(data_names_array(speed_index -1), 'Location','southwest');

figure;
plot(data_array(height_index, :), data_array(alpha_index, :), 'b');
xlabel('Height (m)');
ylabel('Angle of Attack (Degrees)');
grid on;
xlim([-2000, max(data_array(height_index, :)) + 2000]);
ylim([5, max(data_array(alpha_index, :)) + 5]);
set(gca, 'xdir', 'reverse');
title('Angle of Attack of a HRV (3DOF)');
legend(data_names_array(height_index -1), 'Location','southwest');
