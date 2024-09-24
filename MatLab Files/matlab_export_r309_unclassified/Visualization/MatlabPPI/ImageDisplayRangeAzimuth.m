% Initial Class for Range (y) Azimuth (x) Display
classdef ImageDisplayRangeAzimuth < ImageDisplayBase
  properties
    m_hSweepPatch;
    m_sweepGradientAngle = 3*pi/180;
  end
  methods
    function obj = ImageDisplayRangeAzimuth(nPixelsY,nPixelsX,maxRange)
      % call base class
      obj = obj@ImageDisplayBase(nPixelsY,nPixelsX,0,2*pi,0,maxRange);
      % assign properties
      % set labels
      set(get(obj.m_hAxisHandle,'XLabel'),'String','Scan Azimuth (deg) ');
      set(get(obj.m_hAxisHandle,'YLabel'),'String','Slant Range (km)');
      % set default axis
      obj.updateAzimuthTicks(45);
      obj.updateRangeTicks(round(maxRange/5));
      % sweep line
      obj.m_hSweepPatch = patch();
      set(obj.m_hSweepPatch, ...
        'XData',[ 0 0 0 0], ...
        'YData',[ 0 0 0 0], ...
        'Clipping','off', ...
        'UserData', maxRange, ...
        'FaceColor',[0 1 0]);
    end
    function updateAzimuthTicks(this,degStep)
      this.updateHorizontalTicks(pi/180*degStep,180/pi);
    end
    function updateRangeTicks(this,rangeStep)
      this.updateVerticalTicks(rangeStep,1e-3);
    end
    function setPixel(this,azRad,rngM,RGB) % intensity 0-255
      azRad = mod(azRad + 4*pi,2*pi); % ensure angle between 0 and 2pi
      this.setPixel@ImageDisplayBase(azRad,rngM,RGB);
    end
    function setSweepPosition(this,azRad)
      azRad = mod(azRad + 4*pi,2*pi);
      x1 = azRad - this.m_sweepGradientAngle;
      maxr = get(this.m_hSweepPatch,'UserData');
      set(this.m_hSweepPatch, ...
        'XData',[x1 x1 azRad azRad],...
        'YData',[0 maxr maxr 0]);
    end
    function setPatchOverlay(this,ID,az,range,azhalfspan,rangehalfspan,RGBONE,thickness) %input in eu
      this.setPatchOverlay@ImageDisplayBase(ID, ...
        az-azhalfspan, ...
        az+azhalfspan, ...
        range-rangehalfspan, ...
        range+rangehalfspan, ...
        RGBONE,thickness);
    end
  end
end

