classdef ImageDisplayBase < handle
  properties
    m_nPixelsX;
    m_nPixelsY;
    m_hFigureHandle;
    m_hAxisHandle;
    m_hImageHandle;
    m_axisLineBrightness = 0.4;
    m_max_x_row_engineeringUnit;
    m_max_y_col_engineeringUnit;
    m_min_x_row_engineeringUnit;
    m_min_y_col_engineeringUnit;
    m_x_row_eu_per_pixel;
    m_y_col_eu_per_pixel;
    m_x_row_vert_tick_scale_factor;
    m_y_col_horz_tick_scale_factor;
    m_imageData;
    m_phosphorGreenSecondsToClear;
    m_phosphorRedSecondsToClear;
    m_phosphorBlueSecondsToClear;
    m_phosphorFrameSubtract;
    m_updateFramesPerSecond;
    m_isReady;
    m_activeOverlayPatches;
    m_activePlotObjects;
    m_activeTextObjects;
    m_masterDrawingTimer = timer( ...
      'Name','MasterTimer',   ...
      'ExecutionMode','fixedRate', ...
      'BusyMode','queue',...
      'Period',0.1 ...
      );
    m_aviFile;

  end
  methods
    function obj = ImageDisplayBase(nPixelsY,nPixelsX,minUnitY,maxUnitY,minUnitX,maxUnitX)
      
      % ass ign properties
      obj.m_phosphorGreenSecondsToClear  = 3;
      obj.m_phosphorRedSecondsToClear = 3;
      obj.m_phosphorBlueSecondsToClear = 3;
      obj.m_updateFramesPerSecond = 10;
      obj.m_nPixelsX = nPixelsX;
      obj.m_nPixelsY = nPixelsY;
      obj.m_max_x_row_engineeringUnit = maxUnitX;
      obj.m_max_y_col_engineeringUnit = maxUnitY;
      obj.m_min_x_row_engineeringUnit = minUnitX;
      obj.m_min_y_col_engineeringUnit = minUnitY;
      obj.m_x_row_eu_per_pixel = (maxUnitX-minUnitX) / nPixelsX;
      obj.m_y_col_eu_per_pixel = (maxUnitY-minUnitY) / nPixelsY;

      obj.m_imageData = uint8(zeros(nPixelsX,nPixelsY,3));

      % set green initially (turn on effect)
      obj.m_imageData(:,:,2) = 128;
      obj.calculatePhosphorFrameSubract();
      % create figure
      obj.m_hFigureHandle = figure;
      set(obj.m_hFigureHandle,'DoubleBuffer','on', ...
        'Renderer','painters');
      set(obj.m_hFigureHandle,'Color',[0 0 0]);
      % display image
      obj.m_hImageHandle = image(obj.m_imageData);
      obj.m_hAxisHandle = get(obj.m_hImageHandle,'Parent');
      set(obj.m_hAxisHandle, ...
        'XColor',[0 obj.m_axisLineBrightness 0], ...
        'YColor',[0 obj.m_axisLineBrightness 0], ...
        'XGrid','on',  ...
        'YGrid','on',  ...
        'GridLineStyle','-', ...
        'YDir','normal', ...
        'NextPlot', 'add', ...
        'ALimMode', 'manual', ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'ZLimMode', 'manual' ...
        );

      % set up image coordinates
      set(obj.m_hImageHandle, ...
        'XData',[minUnitY maxUnitY], ...
        'YData',[minUnitX maxUnitX] ...
        );
      axis tight;
     

      
      %'Units','pixels' ...
      % set the timer for decay
      % register w/ timer
      try
        stop(obj.m_masterDrawingTimer);
      catch
        obj.m_masterDrawingTimer = timer( ...
          'Name','MasterTimer', ...
          'ExecutionMode','fixedRate',  ...
          'BusyMode','queue',...
          'Period',0.1 ...
          );
      end
      if ( iscell(obj.m_masterDrawingTimer.UserData) )
        temp = obj.m_masterDrawingTimer.UserData;
        set(obj.m_masterDrawingTimer,   ...
          'TimerFcn',@obj.masterTimerFunction, ...
          'ErrorFcn',@obj.timerError,  ...
          'Name','MasterDrawingTimer', ...
          'UserData', {temp{:},obj});
      else
        set(obj.m_masterDrawingTimer, ...
          'TimerFcn',@obj.masterTimerFunction, ...
          'ErrorFcn',@obj.timerError, ...
          'Name','MasterDrawingTimer', ...
          'UserData', {obj});
      end
      start(obj.m_masterDrawingTimer);
      % kick off the timer
      % add callback to figure that manages the timer if the window is
      % closed.
      set(obj.m_hFigureHandle,'DeleteFcn',@obj.figureCloseCleanup);
      this.m_aviFile = [];
      obj.m_isReady = 1;
    end
    
     function masterTimerFunction(varargin)
       displayHandles = varargin{1}.m_masterDrawingTimer.UserData;
       for nDisplay = 1:length(displayHandles)
         displayHandles{nDisplay}.updateImageData();
       end
       drawnow;
     end    
    
    function updateImageData(this)
      this.m_imageData(:,:,1) = this.m_imageData(:,:,1)-this.m_phosphorFrameSubtract(1);
      this.m_imageData(:,:,2) = this.m_imageData(:,:,2)-this.m_phosphorFrameSubtract(2);
      this.m_imageData(:,:,3) = this.m_imageData(:,:,3)-this.m_phosphorFrameSubtract(3);
      set(this.m_hImageHandle,'CData',this.m_imageData);
    end
    function updateHorizontalTicks(this,horizStepEU,horizScaleFactor)
      if ( nargin == 2 )
        horizScaleFactor = 1.0;
      end
      set(this.m_hAxisHandle, ...
        'XTick',this.m_min_y_col_engineeringUnit:horizStepEU:this.m_max_y_col_engineeringUnit);
      set(this.m_hAxisHandle, ...
        'XTickLabel',horizScaleFactor*get(this.m_hAxisHandle,'XTick'));
    end
    function updateVerticalTicks(this,verticalStepEU,verticalScaleFactor)
      if ( nargin == 2 )
        verticalScaleFactor = 1.0;
      end
      set(this.m_hAxisHandle, ...
        'YTick',this.m_min_x_row_engineeringUnit:verticalStepEU:this.m_max_x_row_engineeringUnit);
      set (this.m_hAxisHandle, ...
        'YTickLabel',verticalScaleFactor*get(this.m_hAxisHandle, 'YTick'));
    end

    function timerError(varargin)
      try
        stop(varargin{1}.m_masterDrawingTimer);
        delete(varargin{1}.m_masterDrawingTimer);
      catch
      end
      varargin{1}.m_isReady = 0;
    end
    function setTimeToDecayRed(this,timeSec)
      this.m_phosphorRedSecondsToClear = timeSec;
      this.calculatePhosphorFrameSubract();
    end
    function setTimeToDecayGreen(this,timeSec)
      this.m_phosphorGreenSecondsToClear = timeSec;
      this.calculatePhosphorFrameSubract();
    end
    function setTimeToDecayBlue(this,timeSec)
      this.m_phosphorBlueSecondsToClear = timeSec;
      this.calculatePhosphorFrameSubract();
    end
    function setGraphicsUpdateFramesPerSec(this,nFrames)
      try
        stop(this.m_masterDrawingTimer);
      catch
        disp('Stale timer object found! Creating a new one.');
        %             this.m_phosphorTimer = timer( ...
        %             'TimerFcn',@this.decayData, ...
        %             'ErrorFcn',@this.timerError, ...
        %             'Name','PhosphorTimer', ...
        %             'ExecutionMode','fixedRate', ...
        %             'Period',round(1000/this.m_updateFramesPerSecond)/1000 ...
        %             );
      end
      this.m_updateFramesPerSecond = nFrames;
      this.calculatePhosphorFrameSubract();
      set(this.m_masterDrawingTimer, ...
        'Period',round(1000/this.m_updateFramesPerSecond)/1000);
      start(this.m_masterDrawingTimer);
    end
    function calculatePhosphorFrameSubract(this)
      this.m_phosphorFrameSubtract = [ ...
        255/(this.m_updateFramesPerSecond*this.m_phosphorRedSecondsToClear) ...
        255/(this.m_updateFramesPerSecond*this.m_phosphorGreenSecondsToClear) ...
        255/(this.m_updateFramesPerSecond*this.m_phosphorBlueSecondsToClear) ];
    end
    function figureCloseCleanup(varargin)
      varargin{1}.m_hFigureHandle = [];
      varargin{1}.m_isReady = 0;
      delete(varargin{1});
    end

    function pokeGreen(this)
      this.m_imageData(:,:,2) = 255*rand(this.m_nPixelsX,this.m_nPixelsY);
    end
    function setPixel(this,y_horiz_col_value_eu, x_vert_row_value_eu,RGB) % intensity 0-255
      % convert range to pixels
      pixel_column = round((y_horiz_col_value_eu-this.m_min_y_col_engineeringUnit) / this.m_y_col_eu_per_pixel);
      % convert elevation to pixels
      pixel_row = round((x_vert_row_value_eu-this.m_min_x_row_engineeringUnit) / this.m_x_row_eu_per_pixel);
      % lower bounds checking
      if(pixel_column == 0)
        pixel_column = 1;
      end
      if(pixel_row == 0)
        pixel_row = 1;
      end
      % upper bounds checking
      if(pixel_column > this.m_nPixelsY)
        pixel_column = this.m_nPixelsY;
      end
      if(pixel_row > this.m_nPixelsX)
        pixel_row = this.m_nPixelsX;
      end
      % set pixel
      this.m_imageData(pixel_row,pixel_column,:) = RGB;
    end
    function setPatchOverlay(this,ID,horzstart,horzstop,vertstart,vertstop,RGBONE,thickness) %input in eu
      % check if a patch of this ID is already open
      hExisting = findobj(this.m_activeOverlayPatches,'flat','Tag',ID);
      %        found = 0;
      %        for k = 1:length(this.m_activeOverlayPatches)
      %             if(st rcmp(get(this.mactiveOve rlayPatches(k),'Tag'),ID))
      %                 found = k;
      %                 break;
      %             end
      %        end
      if ( isempty(hExisting) )
        hP = patch( ...
          [horzstart horzstart horzstop horzstop], ...
          [vertstart vertstop vertstop vertstart], ...
          'w',   ...
          'FaceColor','none', ...
          'Clipping', 'off', ...
          'Tag',ID, ...
          'Parent', this.m_hAxisHandle, ...
          'LineWidth', thickness, ...
          'EdgeColor', RGBONE);
        this.m_activeOverlayPatches = [ this.m_activeOverlayPatches hP ];
      else

        set(hExisting,   ...
          'EdgeColor',RGBONE, ...
          'XData', [horzstart horzstart horzstop horzstop], ...
          'YData', [vertstart vertstop vertstop vertstart] ...
          );
      end
    end

    function removePatchOverlay(this,ID)
      hExisting = findobj(this.m_activeOverlayPatches,'Tag',ID);
      if ( ~isempty(hExisting) )
        j = find(this.m_activeOverlayPatches == hExisting);
        if ( j == 1 )
          this.m_activeOverlayPatches = this.m_activeOverlayPatches(2:end);
        elseif ( j == length(this.m_activeOverlayPatches) )
          this.m_activeOverlayPatches = this.m_activeOverlayPatches(1:(end-1));
        else
          this.m_activeOverlayPatches = this.m_activeOverlayPatches([1:(j-1) (j+1):end]);
        end
        delete(hExisting);
      end
    end
    function setPlotOverlay(this,ID,horizontal,vertical,RGBONE,markerstring,markersize)
      hExisting = findobj(this.m_activePlotObjects,'Tag',ID);
      if ( isempty(hExisting) )
        hP = plot( ...
          this.m_hAxisHandle, ...
          horizontal,vertical, ...
          'Marker',markerstring, ...
          'MarkerSize',markersize, ...
          'Clipping','off', ...
          'Tag', ID, ...
          'Color',RGBONE);
        this.mactivePlotObjects = [ this.m_activePlotObjects hP ];
      else
        set(hExisting, ...
          'Marker', markerstring, ...
          'MarkerSize', markersize, ...
          'Color', RGBONE, ...
          'XData', horizontal, ...
          'YData', vertical ...
          )
      end
    end
    function removePlotOverlay(this,ID)
      hExisting = findobj(this.mactivePlotObjects,'Tag',ID);
      if ( ~isempty(hExisting) )
        j = find(this.m_activePlotObjects == hExisting);
        if ( j == 1)
          this.m_activePlotObjects = this.m_activePlot0bjects(2:end);


        elseif ( j == length (this.m_activePlotObjects) )
          this.m_activePlotObjects = this.m_activePlotObjects(1: (end-1));
        else
          this.m_activePlotObjects = this.m_activePlotObjects([1:(j-1) (j+1):end]);
        end
        delete(hExisting);
      end
    end
    function setTextObject(this,ID,horizontal,vertical,string,RGBONE,size)
      hExisting = findobj(this.m_activeTextObjects,'Tag',ID);
      if ( isempty(hExisting) )
        hT = text(horizontal,vertical,string, ...
          'Parent',this.m_hAxisHandle, ...
          'FontSize',size,...
          'Tag', ID, ...
          'Color',RGBONE);
        this.m_activeTextObjects = [ this.m_activeTextObjects hT ];
      else
        set(hExisting, ...
          'Color', RGBONE, ...
          'FontSize',size, ...
          'Position',[horizontal vertical 0] ...
          );
      end
    end
    function removeTextObject(this,ID)
      hExisting = findobj(this.m_activeTextObjects,'Tag',ID);
      if ( ~isempty(hExisting) )
        j = find(this.m_activeTextObjects == hExisting);
        if ( j == 1 )
          this.m_activeTextObjects = this.m_activeTextObjects(2:end);
        elseif ( j == length(this.m_activeTextObjects) )
          this.mactiveTextObjects = this.m_activeTextObjects(1:(end-1));
        else
          this.m_activeTextObjects = this.m_activeTextObjects([1:(j-1) (j+1):end]);
        end
        delete(hExisting);
      end
    end
    function delete(this)
      try
        stop(this.m_masterDrawingTimer);
        delete(this.m_masterDrawingTimer);
      catch
        disp('Could not stop or delete the master drawing timer! ');
      end
      try
        if(~isempty(this.m_hFigureHandle))
          close(this.m_hFigureHandle);
        end
      catch
        disp('Could not close figure window!');
      end
    end
    
    function setAviFile(this,stringName,fps,compressionString,quality,keyframeTime)
        this.m_aviFile                  = avifile(stringName,'Compression',compressionString);
        this.m_aviFile.fps              = fps;
        this.m_aviFile.Quality          = quality;
        this.m_aviFile.KeyFramePerSec   = keyframeTime;
    end
    
    function setAviRecordOn(this,varargin)
        % stop the timer
        stop(this.m_masterDrawingTimer);
        if(nargin > 1)
            aviFramesToRecord = varargin{1};
        else
            aviFramesToRecord = 1e16;
            disp('Ctrl-C to stop');
        end
        for k = 1:aviFramesToRecord
            disp(['Making Frame ' num2str(k)]);
            if(nargin > 2)
                varargin{2}(k);
            end
            this.masterTimerFunction(this);
            this.m_aviFile = addframe(this.m_aviFile,this.m_hFigureHandle);
        end
        this.m_aviFile = close(this.m_aviFile);
        start(this.m_masterDrawingTimer);
    end
  end
end
