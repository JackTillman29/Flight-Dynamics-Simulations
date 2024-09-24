% Initial Class for Elevation (y) Range (x) Display
classdef ImageDisplayElevationRange < ImageDisplayBase
    properties
    end
    methods
        function obj = ImageDisplayElevationRange(nPixelsX,nPixelsY,maxRange,maxElevationDeg)
            % Call base class
            obj = obj@ImageDisplayBase(nPixelsY,nPixelsX,0,maxRange,0,pi/180*maxElevationDeg);
            % assign properties
            % set labels
            set(get(obj.m_hAxisHandle,'YLabel'),'String','Scan Elevation (deg)');
            set(get(obj.m_hAxisHandle,'XLabel'),'String','Slant Range (km)');
            % set default axis
            obj.updateElevationTicks(10);
            obj.updateRangeTicks(round(maxRange/5));
        end
        function updateElevationTicks(this,degStep)
            this.updateVerticalTicks(pi/180*degStep,180/pi);
        end
        function updateRangeTicks(this,rangeStep)
            this.updateHorizontalTicks(rangeStep,1e-3);
        end
        
        function setPixel(this,rangeMeters,elevationRad,RGB) % intensity 0-255
            this.setPixel@ImageDisplayBase(rangeMeters,elevationRad,RGB);
        end

        function setPatchOverlay(this,ID,range,el,rangehalfspan,elhalfspan,RGBONE,thickness) % input in eu
            this.setPatchOverlay@ImageDisplayBase(ID, ...
                range-rangehalfspan, ...
                range+rangehalfspan, ...
                el-elhalfspan, ...
                el+elhalfspan, ...
                RGBONE,thickness);
        end
    end
end
