classdef CartesianFrame < handle
    % Transformation matrix properties:
    % This frame definition in basis vector set is
    % [ --- ux --- ]
    % [ --- uy --- ]
    % [ --- uz --- ]
    %
    % Basis vector set in this frame is
    % [ |  |  |  ]
    % [ ux uy uz ]
    % [ |  |  |  ]
    properties
        ux = [1 0 0]';
        uy = [0 1 0]';
        uz = [0 0 1]';
        name = 'UNK';
    end
    methods
        function obj = CartesianFrame(name,defn)
            if( ~exist('defn','var'))
                defn = eye(3);
            end
            obj.ux = defn(1,:);
            obj.uy = defn(2,:);
            obj.uz = defn(3,:);
            obj.name = name;
        end
        function RelativeToFrame(objIn,M)
            objIn.ux = M * objIn.ux;
            objIn.uy = M * objIn.uy;
            objIn.uz = M * objIn.uz;
        end
        function R = AsMatrix(objIn)
            R = [objIn.ux; objIn.uy; objIn.uz];
        end
        function draw(objIn,origin,scaleFactor)
            if(~exist('origin'))
                origin = [0 0 0];
            end
            if(~exist('scaleFactor'))
                scaleFactor = 1;
            end
            hLines = [ ...
                line(scaleFactor*[0 objIn.ux(1)]+origin(1),scaleFactor*[0 objIn.ux(2)]+origin(2),scaleFactor*[0 objIn.ux(3)]+origin(3)) ...
                line(scaleFactor*[0 objIn.uy(1)]+origin(1),scaleFactor*[0 objIn.uy(2)]+origin(2),scaleFactor*[0 objIn.uy(3)]+origin(3)) ...
                line(scaleFactor*[0 objIn.uz(1)]+origin(1),scaleFactor*[0 objIn.uz(2)]+origin(2),scaleFactor*[0 objIn.uz(3)]+origin(3))];
            set(hLines(1),'Color','r','LineWidth',2);
            set(hLines(2),'Color','g','LineWidth',2);
            set(hLines(3),'Color','b','LineWidth',2);
            ht(1)=text(scaleFactor*objIn.ux(1),scaleFactor*objIn.ux(2)+origin(2),scaleFactor*objIn.ux(3)+origin(3),[objIn.name ' X']);
            ht(2)=text(scaleFactor*objIn.uy(1)+origin(1),scaleFactor*objIn.uy(2)+origin(2),scaleFactor*objIn.uy(3)+origin(3),[objIn.name ' Y']);
            ht(3)=text(scaleFactor*objIn.uz(1)+origin(1),scaleFactor*objIn.uz(2)+origin(2),scaleFactor*objIn.uz(3)+origin(3),[objIn.name ' Z']);
            set(ht,'FontSize',6);
            
            % draw planes
            da = 2*pi / 30;
            ang = 0:da:(2*pi);
            xpts = cos(ang);
            ypts = sin(ang);
            
            M = AsMatrix(objIn);
            temp = M' * scaleFactor*[xpts;ypts;0*xpts;];

            hx=patch(temp(1,:)+origin(1),temp(2,:)+origin(2),temp(3,:)+origin(3),'b');
            
            temp = M' * scaleFactor*[xpts;0*ypts;ypts;];
            hy=patch(temp(1,:)+origin(1),temp(2,:)+origin(2),temp(3,:)+origin(3),'g');
            
            temp = M' * scaleFactor*[0*xpts;xpts;ypts;];
            hz=patch(temp(1,:)+origin(1),temp(2,:)+origin(2),temp(3,:)+origin(3),'r');
            
            
            set([hx hy hz],'EdgeAlpha',0.2,'FaceAlpha',0.1);
            
        end
        function disp(objIn)
            fprintf('|%+5.3f %+5.3f %+5.3f|\n|%+5.3f %+5.3f %+5.3f|\n|%+5.3f %+5.3f %+5.3f|\n|', ...
                objIn.ux(1),objIn.ux(2),objIn.ux(3), ...
                objIn.uy(1),objIn.uy(2),objIn.uy(3), ...
                objIn.uz(1),objIn.uz(2),objIn.uz(3) );
                
        end
    end
end