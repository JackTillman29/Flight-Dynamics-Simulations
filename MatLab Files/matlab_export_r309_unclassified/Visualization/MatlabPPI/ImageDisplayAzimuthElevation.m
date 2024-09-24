%Initial Class for Range (y) Azimuth (x) Display
classdef ImageDisplayAzimuthElevation < ImageDisplayBase
    properties
        m_hSweepPatch;
        m_sweepGradientAngle = 3*pi/180;
    end
    methods
        function obj = ImageDisplayAzimuthElevation(nPixelsY,nPixelsX,maxElevation)
            % call base class
             obj = obj@ImageDisplayBase(nPixelsY,nPixelsX,0,2*pi,0,pi/180*maxElevation);
            % assign properties
            % set labels
             set(get(obj.m_hAxisHandle,'XLabel'),'String','Scan Azimuth (deg)');
             set(get(obj.m_hAxisHandle,'YLabel'),'String','Elevation (deg)');
            % set default axis
            obj.updateAzimuthTicks(45);
            obj.updateElevationTicks(10);
            % sweep line
            obj.m_hSweepPatch = patch();
            set(obj.m_hSweepPatch, ...
                 'XData', [ 0 0 0 0], ...
                 'YData', [ 0 0 0 0], ...
                 'Clipping','off', ...
                 'UserData',pi/180*maxElevation, ...
                 'FaceColor',[0 1 0]);
        end
        function updateAzimuthTicks(this,degStep)
            this.updateHorizontalTicks(pi/180*degStep,180/pi);
        end
        function updateElevationTicks(this,degStep)
            this.updateVerticalTicks(pi/180*degStep,180/pi);
        end
        function setPixel(this,azRad,elRad,RGB) % intensity 0-255
            azRad = mod(azRad + 4*pi,2*pi);
            this.setPixel@ImageDisplayBase(azRad,elRad,RGB);
        end
        function setSweepPosition(this,azRad)
            azRad = mod(azRad + 4*pi,2*pi);
            x1 = azRad - this.m_sweepGradientAngle;
            elMax = get(this.m_hSweepPatch,'UserData');
            set(this.m_hSweepPatch, ...
                 'XData',[x1 x1 azRad azRad],...
                 'YData',[0 elMax elMax 0]);
        end
        function setPatchOverlay(this, ID, az, el, azhalfspan, elhalfspan, RGBONE, thickness) % input in eu
             this.setPatchOverlay@ImageDisplayBase(ID, ...
                 az-azhalfspan, ...
                 az+azhalfspan, ...
                 el-elhalfspan, ...
                 el+elhalfspan, ...
                 RGBONE,thickness);                 
        end
     end
end