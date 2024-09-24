close all; clear all; clc

% RED DOT and BLUE CIRCLE SHOULD BE AT THE SAME PLACE
%  LINE ANNOTATION SHOULD POINT FROM BOTTOM-LEFT TO THE RED DOT/BLUE CIRCLE

N = 1000
x = 1:N;
y = randn(1,N);

hfig = figure;

for k = 300:20:1000
plot(x,y);
hold on
hax = gca;

idx = find(x == k);

pt_in_axes = [k y(idx)];

pt_in_fig = pointInAxesToFigurePosition(pt_in_axes,hax);

plot(pt_in_axes(1),pt_in_axes(2),'r.','MarkerSize',15)

hannot = annotation('line',[0 pt_in_fig(1)],[0 pt_in_fig(2)]);


pt_in_axes2 = pointInFigureToAxesPosition(pt_in_fig,hax);

plot(pt_in_axes2(1),pt_in_axes2(2),'bo','MarkerSize',15)
hold off

pause(0.1)

delete(hannot)

end



%%

