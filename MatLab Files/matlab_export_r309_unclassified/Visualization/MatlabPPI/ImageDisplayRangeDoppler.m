%Initial Class for Range (y) Azimuth (x) Display
classdef ImageDisplayRangeDoppler < ImageDisplayBase
    properties
        m_nRangeBins;
        m_nDopplerFilters;
    end
    methods
        function obj = ImageDisplayRangeDoppler(nRangeBins,nDopplerFilters)
            % call base class
             obj = obj@ImageDisplayBase(nDopplerFilters,nRangeBins,1,nDopplerFilters,1,nRangeBins);
            % assign properties
            % set labels
             set(get(obj.m_hAxisHandle,'XLabel'),'String','Doppler Filter (#)');
             set(get(obj.m_hAxisHandle,'YLabel'),'String','Range Bin (#)');
            % set default axis
            obj.updateRangeTicks(nRangeBins);
            obj.updateDopplerTicks(nDopplerFilters);
            set(obj.m_hAxisHandle,'GridLineStyle','none');
        end
        function updateDopplerTicks(this,degStep)
            this.updateHorizontalTicks(1,1);
        end
        function updateRangeTicks(this,degStep)
            this.updateVerticalTicks(1,1);
        end
        function setPatchOverlay(this,dopplerFilter, rangeBin, dopplerWidth, rangeWidth, RGBONE, thickness) % input in eu
             this.setPatchOverlay@ImageDisplayBase(ID, ...
                 dopplerFilter-dopplerWidth, ...
                 dopplerFilter+dopplerWidth, ...
                 rangeBin-rangeWidth, ...
                 rangeBin+rangeWidth, ...
                 RGBONE,thickness);                 
        end
     end
end