classdef vidSyncObj < handle
    properties
        vidObj;
        vidReadObj;
        imAxes;
        parentFigure;
        imHandle;
        frameRate;
        filename;
        nFrames;
        framesRemaining;
        nMasterFrames = [];
        drawDecimateFactor = 1;
        drawDecimateCounter = 0;
        
    end
    methods
         function this = vidSyncObj(filename,parentFig,masterFps)
             this.filename = filename;
             this.vidReadObj = VideoReader(filename);
             this.nFrames = this.vidReadObj.FrameRate * this.vidReadObj.Duration;
             this.frameRate = this.vidReadObj.FrameRate;
             this.framesRemaining = this.nFrames;
             this.imHandle = [];
             this.parentFigure = parentFig;
             this.imAxes = [];
             
             disp(['Frame rate factor is ' num2str(masterFps/this.frameRate)]);

             
         end
         
         function this = setup(this,position,drawDecimate)
             figure(this.parentFigure);
             this.imAxes = axes();
             draw_data = this.vidReadObj.readFrame();
             this.framesRemaining = this.framesRemaining - 1;
             this.imHandle = image(draw_data,'Parent',this.imAxes);
             set(this.imAxes,'Visible','off', ...
                 'DataAspectRatio',[1 1 1]);
             set(this.imAxes,'Position',position);
             this.drawDecimateFactor = drawDecimate;
             
             this.nMasterFrames = this.drawDecimateFactor * this.nFrames;
         end
         
         function this = draw(this)
             % first increment decimate counter
             this.drawDecimateCounter = this.drawDecimateCounter + 1;
             
             % now check if it is time
             if(this.drawDecimateCounter == this.drawDecimateFactor)
                 if(this.framesRemaining > 0)
                     draw_data = this.vidReadObj.readFrame();
                     set(this.imHandle,'CData',draw_data);
                     this.framesRemaining = this.framesRemaining - 1;
                     this.drawDecimateCounter = 0;
                 else
                     disp(['End of file: ' this.filename]);
                 end
             end
         end
    end
end