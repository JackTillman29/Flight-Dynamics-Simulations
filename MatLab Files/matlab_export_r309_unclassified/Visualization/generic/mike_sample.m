% say you have a beam that is 3degx5deg, pointed at az=45 and el = 10.
% you want the beam to be blue, and see-through

% first, generate the xy data for the beam
axis equal;
xlim([0 180]);
ylim([0 90]);
[beamx,beamy]=ellipse_patch( ...
    45, ...     % farfield az 
    10, ...     % farfield el
    3, ...      % az beamwidth
    5, ...      % el beamwidth
    50);        % # of points in the ellipse (aesthetics)

h(1) = patch(beamx,beamy,'k');
set(h(1),'FaceColor','b','FaceAlpha',0.4,'EdgeAlpha',0.1);

pause(0.5);

[beamx,beamy]=ellipse_patch( ...
    47, ...     % farfield az 
    10, ...     % farfield el
    3, ...      % az beamwidth
    5, ...      % el beamwidth
    50);        % # of points in the ellipse (aesthetics)

h(2) = patch(beamx,beamy,'k');
set(h(2),'FaceColor','r','FaceAlpha',0.4,'EdgeAlpha',0.1);

pause(0.5);

[beamx,beamy]=ellipse_patch( ...
    49, ...     % farfield az 
    10, ...     % farfield el
    3, ...      % az beamwidth
    5, ...      % el beamwidth
    50);        % # of points in the ellipse (aesthetics)

h(3) = patch(beamx,beamy,'k');
set(h(3),'FaceColor','g','FaceAlpha',0.4,'EdgeAlpha',0.1);

pause(0.5);

delete(h(1));

pause(0.5);

delete(h(2));

pause(0.5);

delete(h(3));


