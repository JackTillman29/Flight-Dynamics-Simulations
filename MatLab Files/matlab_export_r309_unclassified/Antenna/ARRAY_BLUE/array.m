classdef array < handle
    properties
        ENUM_MODE_UNIFORM = 0;
        ENUM_MODE_AMPFILE = 1;  %associated with rx amp file
        ENUM_MODE_DAZ_AMPFILE = 2;
        ENUM_MODE_DEL_AMPFILE = 3;
        ENUM_MODE_ABF = 4;
        ENUM_MODE_TXAMPFILE = 5; %associated with tx amp file, if specified
        pos
        elem_area
        elemfac_gain
        amp,ampel,ampaz,amptx
        element_off
        rng
        rng_feed
        phs
        phs_steer
        phs_sphcmp
        abf_term
        dly
        wavelength
        efficiency
        cos_pow_scn = 1.2
        twopi_ovr_lambda
        neScale
        nElements,elemFacPresent
        vPeak,aVar,pVar,azSteer,elSteer,azSteerThresh,elSteerThresh
        adcSample
        pShiftQuantization,needPhaseQuantized
        pMin,pMax,pBits,pStep
        isInit = 0 %as with Fortran, set flag to only allow init once
        upd_weightsScaled = 0; %as with Fortran, set flag to only allow scaling once
        
        elementPattern %the normalized radiation pattern of one element from the array
        elementPatternCross
        elementPolCart = struct('U',[],'V',[],'W',[]); %Contains polarzation of radiation pattern in cartiesan parts
        elemAz, elemEl
        
        acBodyToAntenna = eye(3); % assume antenna straight out nose
        
        % matlab specific
        calcMethod = 2; % 2 = vectorized, 1 = exact FORTRAN copy
        gpuMethod = 0;
        precisionClass = 'single';
        name = '';
        plot_marker_size = 32;
        plot_2d_res_deg = 1.0;
        
        upd_isvalid = 0;
        upd_elementPattern = [];
        upd_elementPatternAzG = [];
        upd_elementPatternElG = [];
        upd_elementPatternPresent = 0;
        upd_elementGaindBV = [];
        upd_ok2updateNoise = 1;
        upd_wn = 1; % array elemental noise vector (or 1 if no noise)
        upd_subArrayNumber = []; % Per element [N x 1] sub array assignment
        upd_subArrayMatrix = []; % Linear algebra compliant sub-array (generate via upd_genSubArrayMatrix)
        %upd_subArrayWgt = [];    % Apply sub-array weighting to full array
        upd_subArrayPos = [];
        % These variables are just for tool robustness and speed of
        % plotting based on user dynamic updates.
        upd_farField_staticArray = [];  % used to store Az/El phasing static portion of array calculation
        upd_farField_staticElement = []; % used to store Az/El phasing for static element data
        upd_farField_staticAzEl = []; % holds static nAz,nEl
        upd_wgt_name = [];
        upd_handle_2d_plots = [-1 -1];

    end
    properties (SetAccess=private)
        upd_subArrayWgt = [];
    end
    methods
        function obj = array(name)
            % empty constructor
            if(~exist('name'))
                obj.name = '';
            else
                obj.name = name;
            end
            
            % check for GPU
            try
                nGPUs = gpuDeviceCount();
            catch
                warning('Could not determine GPU count. Defaulting to CPU based computation.');
                nGPUs = 0;
            end
            if(nGPUs > 0)
                fprintf('Found GPU!\n ** Call <array>.useGPU() to use GPU\n');
            end
            obj.aVar = 0.0;
            obj.pVar = 0.0;
        end
        
        %         function setElementFactor(this,angleRad,gainDbi)
        %             this.elemfac_gain = 10.^(gainDbi/20);
        %             this.elemfac_ang = angleRad;
        %             this.elemFacPresent = 1;
        %         end
        
        function setElementFactor(this,azAng,elAng,gain_v)
            this.elemfac_gain = gain_v;
            this.elemfac_az = azAng;
            this.elemfac_el = elAng;
            this.elemFacPresent = 1;
        end
        
        function turnOff(this,idxs)
            % turnOff(this,idxs), idxs = element #s or 'none'
            this.element_off = 0*this.amp;
            if(ischar(idxs))
                if(strcmpi(idxs,'all'))
                    this.element_off = 0*this.amp + 1;
                elseif(strcmpi(idxs,'none'))
                    this.element_off = 0*this.amp;
                else
                    error(['turnOff: unknown option: ' idxs]);
                end
            else
                this.element_off(idxs) = 1;
            end
            this.neScale = 1.0 / sum(this.element_off == 0);
            try
                this.vPeak = (sqrt( 4 * pi * (this.elem_area/this.wavelength^2) * this.efficiency / this.neScale));
            catch
                warning('Unable to compute vPeak.');
            end
        end
        
        function turnOn(this,idxs)
            % turnOff(this,idxs), idxs = element #s or 'none'
            this.element_off = 0*this.amp+1;
            if(ischar(idxs))
                if(strcmpi(idxs,'all'))
                    this.element_off = 0*this.amp;
                elseif(strcmpi(idxs,'none'))
                    this.element_off = 0*this.amp+1;
                else
                    error(['turnOn: unknown option: ' idxs]);
                end
            else
                this.element_off(idxs) = 0;
            end
            this.neScale = 1.0 / sum(this.element_off == 0);
            this.vPeak = (sqrt( 4 * pi * (this.elem_area/this.wavelength^2) * this.efficiency / this.neScale));
        end
        
        
        function changeInstNormal(this,y,p,r)
            % installation normal (in E,N,U)
            %  + pitch == pointing the array DOWN
            %  - pitch == pointing the array UP
            cy = cos(y);
            cp = cos(p);
            cr = cos(r);
            sy = sin(y);
            sp = sin(p);
            sr = sin(r);
            R = [ ...
                cp*cy               cp*sy               -sp
                sr*sp*cy - cr*sy    cr*cy + sr*sp*sy    sr*cp
                sr*sy + cr*sp*cy    cr*sp*sy - sr*cy    cr*cp ...
                ];
            this.acBodyToAntenna = R;
        end
        function [azAnt,elAnt]=xformAzEl(this,az,el)
            % vectorized + backwards compatibility:
            %      MAKE SURE az AND el ARE ROW VECTORS!
            [ux,uy,uz] = sph2cart(az,el,1.0);
            uant =  this.acBodyToAntenna * [ux;uy;uz];
            if(numel(az) == 1)
                [azAnt,elAnt,temp]=cart2sph(uant(1),uant(2),uant(3));
            else
                [azAnt,elAnt,temp]=cart2sph(uant(1,:),uant(2,:),uant(3,:));
            end
        end
        
        function InitStapArray(this,wavelength,efficiency,posfile,ampfile, ...
                ampfile_az,ampfile_el,elfacfile,txampfile)
            
            if this.isInit 
                disp('[AESA] Array was previously initialized');
                return;
            end
            
            if( ~exist('elfacfile','var') )
                elfacfile = [];
            end
            
            % fortran constructor
            this.pShiftQuantization = 0;
            this.needPhaseQuantized = 0;
            this.wavelength = wavelength;
            this.twopi_ovr_lambda = (2*pi / wavelength);
            this.efficiency = efficiency;
            
            fid = fopen(posfile,'r');
            if(fid == -1)
                error(['Could not open ' posfile]);
            end
            this.nElements = fscanf(fid,'%d',1);
            this.preAllocate();
            
            % Define Array Geometry
            for k = 1 : this.nElements
                this.pos(:,k) = fscanf(fid,'%f',3);
            end
            fclose(fid);
            
            fprintf('=== Array Data ===\nElements: %d\n',this.nElements);
            fprintf('Found %d elements\n',this.nElements);
            
            %recenterGeometry(this);
            
            this.aVar = (0.0);
            this.pVar = (0.0);
            
            this.neScale = (1.0 / this.nElements);
            % jcl(-)this.vPeak = (sqrt(pi * this.efficiency * this.nElements));
            this.vPeak = (sqrt( 4 * pi * (this.elem_area/this.wavelength^2) * this.efficiency * this.nElements));% jcl(+)
            
            if( ~isempty(ampfile) )
                fid = fopen(ampfile,'r');
                nElem = fscanf(fid,'%d',1);
                if nElem ~= this.nElements 
                    error(['inconsistent number of elements in ' ampfile]);
                end
                rawdata = fscanf(fid,'%f');
                fclose(fid);
                if length(rawdata)==nElem                
                    %there is one column of data in the file - real data
                    this.amp = rawdata(:).';
                else
                    %there are two columns of data in the file
                    %first columns is magnitude, second is phase in radians
                   ampmag = rawdata(1:2:end).';
                   ampphs = rawdata(2:2:end).';
                   this.amp = ampmag.*exp(1i*ampphs);
                end
            end
            
            if( ~isempty(ampfile_az) )
                fid = fopen(ampfile_az,'r');
                nElem = fscanf(fid,'%d',1);
                if nElem ~= this.nElements 
                    error(['inconsistent number of elements in ' ampfile_az]);
                end
                rawdata = fscanf(fid,'%f');
                fclose(fid);
                if length(rawdata)==nElem                
                    %there is one column of data in the file - real data
                    this.ampaz = rawdata(:).';
                else
                    %there are two columns of data in the file
                    %first columns is magnitude, second is phase in radians
                   ampmag = rawdata(1:2:end).';
                   ampphs = rawdata(2:2:end).';
                   this.ampaz = ampmag.*exp(1i*ampphs);
                end
            end
            
            
           if( ~isempty(ampfile_el) )
                fid = fopen(ampfile_el,'r');
                nElem = fscanf(fid,'%d',1);
                if nElem ~= this.nElements 
                    error(['inconsistent number of elements in ' ampfile_el]);
                end
                rawdata = fscanf(fid,'%f');
                fclose(fid);
                if length(rawdata)==nElem                
                    %there is one column of data in the file - real data
                    this.ampel = rawdata(:).';
                else
                   %there are two columns of data in the file
                    %first columns is magnitude, second is phase in radians
                   ampmag = rawdata(1:2:end).';
                   ampphs = rawdata(2:2:end).';
                   this.ampel = ampmag.*exp(1i*ampphs);
                end
            end
            
            if( ~isempty(elfacfile) )
                error('Use InitStapArray2DElement for 2D element file');
            end
            
            if nargin == 9  %there is an amplitude file for transmit given
                if( ~isempty(txampfile) )
                    fid = fopen(txampfile,'r');
                    nElem = fscanf(fid,'%d',1);
                    if nElem ~= this.nElements
                        error(['inconsistent number of elements in ' txampfile]);
                    end
                    rawdata = fscanf(fid,'%f');
                    fclose(fid);
                    if length(rawdata)==nElem
                        %there is one column of data in the file - real data
                        this.amptx = rawdata(:).';
                    else
                        %there are two columns of data in the file
                        %first columns is magnitude, second is phase in radians
                        ampmag = rawdata(1:2:end).';
                        ampphs = rawdata(2:2:end).';
                        this.amptx = ampmag.*exp(1i*ampphs);
                    end
                end
            end
                        
            this.rng = sqrt(sum(this.pos.^2,1));
            
            this.isInit = 1;
            
        end % init stap array
        
        function InitStapArray2DElement(this,wavelength,efficiency,...
                posfile,...  %element positions
                ampfile, ...  %receive amplitude weighting
                ampfile_az,....%receive amp weighing diff-az channel
                ampfile_el,...  %receive amp weighing diff-el channel
                elfacfile,...  %element antenna pattern (complex)
                txampfile)     %optional transmit weighting
            
             %assumes efacfile has complex data, with two columns (real
             %left, imaginary right). Second column has zeros for angle information
             %Example:
%             541         <<<--- number of az points
%             361         <<<--- number of el points
%             -1.57079632679489655800  0.00000000000000000000 
%             -1.55334303427495346028  0.00000000000000000000
%             -1.53588974175500991848  0.00000000000000000000
%             etc.
            
            if this.isInit 
                disp('[AESA] Array was previously initialized');
                return;
            end

            if isempty(elfacfile)
                error('No element file specified');
            end

            
            %amplitude weighting on transmit may or may not be present
            if nargin == 8 || (nargin==9 && isempty(txampfile))
                InitStapArray(this,wavelength,efficiency,...
                    posfile,ampfile, ...
                    ampfile_az,ampfile_el,[],[]);
            else
                InitStapArray(this,wavelength,efficiency,...
                    posfile,ampfile, ...
                    ampfile_az,ampfile_el,[],txampfile);
            end
            
            this.upd_elementPatternPresent = 1;
            
            rawdata = dlmread(elfacfile); 
            
            numAz = rawdata(1);
            numEl = rawdata(2);
            
            rawdata_c1 = rawdata(3:end,1); %column 1, angles and real part of data
            rawdata_c2 = rawdata(3:end,2); %column 2, imaginary part of data
            
            this.elemEl = rawdata_c1(1:numEl);
            this.elemAz = rawdata_c1(numEl+1:numEl+1:numEl+numAz+numEl*numAz);
            
            [this.upd_elementPatternAzG,this.upd_elementPatternElG] = ...
                meshgrid(this.elemAz,this.elemEl);
            
            k = 1;
            for j = numEl + 1:numEl+1:numAz+numEl+numAz*numEl
                for kk = 1:numEl
                    %elevation along rows, az along columns
                    this.elementPattern(kk,k) =rawdata_c1(j+kk) + 1i*rawdata_c2(j+kk);
                end
                k = k+1;
            end
            
            this.upd_elementPattern = this.elementPattern;
            
            % Populate peak element gain as scalar value (for reference)
            %if this is unity, this implies a normalized element pattern
            this.upd_elementGaindBV = ...
                10*log10( max(abs(this.upd_elementPattern(:))) );
            
            %normalize the pattern, and apply the peak in upd_getGain
            %this allows for manual control of upd_elementGaindBV after
            %initialization.
             this.upd_elementPattern = ...
                 this.upd_elementPattern/(10^(this.upd_elementGaindBV/10));
            
            if this.upd_elementGaindBV<=1.01
                warning(['Element Pattern is normalized, peak magnitude is: ' ...
                    num2str(this.upd_elementGaindBV) 'dBV']);
                warning('Legacy getGain uses normalized patterns...BUT');
                warning('New upd_getGain expects element gain in the pattern');
                warning('If using upd_getGain, set upd_elementGaindBV manually');
            end
            
        figure;
        surf(180/pi*this.elemAz,180/pi*this.elemEl,20*log10(abs(this.elementPattern)));
        view(2);shading flat;
        colorbar;xlabel('AZ (deg)');ylabel('EL (deg)');
        title({'Element Pattern (dBi)'; ...
            ['Peak Gain (dBi): ' num2str(20*log10(abs(this.upd_elementGaindBV)))]});

        end % init stap array with 2D Element factor file
        
        function useGPU(this)
            if(gpuDeviceCount())
                this.precisionClass = 'gpuArray';
            else
                this.precisionClass = 'single';
                warning('No GPU found. Using CPU single precision.');
            end
        end
        
        function computeFeedPhase(this,focalDistance,feedXOffset,feedYOffset)
            r2 = ...
                (this.pos(1,:) - feedXOffset).^2 + ...
                (this.pos(2,:) - feedYOffset).^2 + ...
                (this.pos(3,:) - focalDistance).^2;
            
            this.rng_feed = sqrt(r2); % distance to the feed el-by-el
            this.phs_sphcmp = this.twopi_ovr_lambda * (this.rng_feed - focalDistance);
            
        end
        
        function setElementXYZ(this,xe,ye,ze,elarea)
            this.nElements = length(xe);
            this.preAllocate();
            this.pos(1,:) = xe;
            this.pos(2,:) = ye;
            this.pos(3,:) = ze;
            %this.recenterGeometry();
            this.neScale = (1.0 / this.nElements);
            this.elem_area = elarea;
        end
        
        function setTxProperties(this,freq,eff)
            this.wavelength = 3e8 / freq;
            this.twopi_ovr_lambda = (2*pi / this.wavelength);
            this.efficiency = eff;
            %this.vPeak = (sqrt(pi * this.efficiency * this.nElements));
            this.vPeak = (sqrt( 4 * pi * (this.elem_area/this.wavelength^2) * this.efficiency * this.nElements));
        end
        
        % This constructor will build a simple square array just for quick
        % testing.
        function CreateRectangularArray(this,lambda,eff,spacefacx,spacefacy,Nx,Ny)
            this.pShiftQuantization = 0;
            this.needPhaseQuantized = 0;
            
            this.twopi_ovr_lambda = (2*pi / lambda);
            this.efficiency = (eff);
            this.wavelength = lambda;
            
            this.nElements = (Nx * Ny);
            xSpace = lambda * spacefacx;
            ySpace = lambda * spacefacy;
            fprintf('=== Array Data ===\nElements: %d\n',this.nElements);
            fprintf('X Element Spacing: %6.4f\nY Element Spacing: %6.4f\n',xSpace,ySpace);
            
            % Define Array Geometry
            this.preAllocate();
            iElement = 1;
            for kx = 1 : Nx
                for ky = 1 : Ny
                    this.pos(:,iElement) = [kx*xSpace ky*ySpace 0]';
                    iElement = iElement + 1;
                end
            end
            %recenterGeometry(this);
            
            this.aVar = (0.0);
            this.pVar = (0.0);
            
            this.neScale = (1.0 / this.nElements);
            %TODO put elemental area in here (vs pi)
            this.elem_area = (lambda*spacefacx) * (lambda*spacefacy);
            %this.vPeak = (sqrt(pi * this.efficiency * this.nElements));
            this.vPeak = (sqrt( 4 * pi * (this.elem_area/this.wavelength^2) * this.efficiency * this.nElements));
            this.turnOn('all');
        end
        
        % This constructor will build a simple square array just for quick
        % testing.
        function CreateArbitraryArray(this,lambda,eff,xSpace,ySpace,zSpace)
            this.pShiftQuantization = 0;
            this.needPhaseQuantized = 0;
            
            this.twopi_ovr_lambda = (2*pi / lambda);
            this.efficiency = (eff);
            this.wavelength = lambda;
            
            this.nElements = length(xSpace);
            %             fprintf('=== Array Data ===\nElements: %d\n',this.nElements);
            %             fprintf('X Element Spacing: %6.4f\nY Element Spacing: %6.4f\nZ Element Spacing: %6.4f\n',...
            %                 xSpace,ySpace,zSpace);
            
            % Define Array Geometry
            this.preAllocate();
            for iElement = 1 : this.nElements
                this.pos(:,iElement) = [...
                    xSpace(iElement) ...
                    ySpace(iElement) ...
                    zSpace(iElement)]';
            end
            %recenterGeometry(this);
            
            this.aVar = (0.0);
            this.pVar = (0.0);
            
            this.neScale = (1.0 / this.nElements);
            %TODO put elemental area in here (vs pi)
            %             this.elem_area = (lambda*spacefacx) * (lambda*spacefacy);
            delta_x = abs(diff(xSpace)); delta_x(delta_x == 0) = inf; minDelta_x = min(delta_x);
            delta_y = abs(diff(ySpace)); delta_y(delta_y == 0) = inf; minDelta_y = min(delta_y);
            %             this.elem_area = 1;
            this.elem_area = minDelta_x * minDelta_y;
            %this.vPeak = (sqrt(pi * this.efficiency * this.nElements));
            this.vPeak = (sqrt( 4 * pi * (this.elem_area/this.wavelength^2) * this.efficiency * this.nElements));
            this.turnOn('all');
        end
        
        function preAllocate(this)
            this.pos = zeros(3,this.nElements,this.precisionClass);
            this.amp = zeros(1,this.nElements,this.precisionClass);
            this.amptx = zeros(1,this.nElements,this.precisionClass);
            this.ampaz = zeros(1,this.nElements,this.precisionClass);
            this.ampel = zeros(1,this.nElements,this.precisionClass);
            this.phs = zeros(1,this.nElements,this.precisionClass);
            this.phs_steer = zeros(1,this.nElements,this.precisionClass);
            this.phs_sphcmp = zeros(1,this.nElements,this.precisionClass);
        end
        
        
        
        
        % This function recenters the element locations such that they
        % are evenly distributed about 0,0,0
        function recenterGeometry(this)
            n = this.nElements;
            
            %xo = 0.5 * (this.pos(1,n) - this.pos(1,1));
            %yo = 0.5 * (this.pos(2,n) - this.pos(2,1));
            %zo = 0.5 * (this.pos(3,n) - this.pos(3,1));
            xo = 0.5 * ( max(this.pos(1,:)) - min(this.pos(1,:)) );
            yo = 0.5 * ( max(this.pos(2,:)) - min(this.pos(2,:)) );
            zo = 0.5 * ( max(this.pos(3,:)) - min(this.pos(3,:)) );
            
            xi = min(this.pos(1,:));
            yi = min(this.pos(2,:));
            zi = min(this.pos(3,:));
            
            for k = 1:n
                this.pos(:,k) = this.pos(:,k) - [xo,yo,zo]';
                this.pos(:,k) = this.pos(:,k) - [xi,yi,zi]';
                % Compute range (assumes the array is centered at 0,0,0)
                %this.rng(k) = DVECTOR_MAGNITUDE(this%pos(:,k))
                
            end
            this.rng = sqrt(sum(this.pos.^2,1));
        end
        
        function wrapPhase(this)
            this.phs = mod(this.phs,2*pi);
            this.phs_steer = mod(this.phs_steer,2*pi);
            this.phs_sphcmp = mod(this.phs_sphcmp,2*pi);
        end
        
        
        function plotGeometry(this,varargin)
            if(this.upd_isvalid == 1)
                xData = this.pos(:,1);
                yData = this.pos(:,2);
                zData = this.pos(:,3);
            else
                xData = this.pos(1,:);
                yData = this.pos(2,:);
                zData = this.pos(3,:);
            end
            % if you pass in a 2nd argument, it will color code
            % the elements.
            %hp = plot3(this.pos(1,:),this.pos(2,:),this.pos(3,:),'.');
            if(nargin == 1)
                hs=scatter3( ...
                    xData, ...
                    yData, ...
                    zData, ...
                    this.plot_marker_size);
                set(hs,'Marker','sq','MarkerFaceColor','r');
            else
                hs=scatter3( ...
                    varargin{1}, ...
                    yData, ...
                    zData, ...
                    this.plot_marker_size, ...
                    varargin{1},'filled');
                set(hs,'Marker','sq');
                colorbar;
            end
            if(this.upd_isvalid == 0)
                view(0,90);
                xlabel('X Direction (m)');
                ylabel('Y Direction (m)');
                zlabel('Z Direction (m)');
                title('Array Geometry');
            else
                set(gca,'YDir','reverse','DataAspectRatio',[1 1 1]);
                ylabel('U (Y), meters');
                zlabel('V (Z), meters');
                xlabel('W (X), Normal, meters');
                title(sprintf('%s Array Geometry',this.name));
                set(hs,'MarkerEdgeColor',0.5*[1 1 1]);
                view(-90,0);
            end
            
        end
        
       function plotElementWeights(this,enum)
            if(this.upd_isvalid == 1)
                xData = this.pos(:,1);
                yData = this.pos(:,2);
                zData = this.pos(:,3);
            else
                xData = this.pos(1,:);
                yData = this.pos(2,:);
                zData = this.pos(3,:);
            end
            % if you pass in a 2nd argument, it will color code
            % the elements.
            %hp = plot3(this.pos(1,:),this.pos(2,:),this.pos(3,:),'.');
            
            if     enum==this.ENUM_MODE_AMPFILE
                wts = this.amp;
            elseif enum == this.ENUM_MODE_DAZ_AMPFILE
                wts = this.ampaz;
            elseif enum == this.ENUM_MODE_DEL_AMPFILE
                wts = this.ampel;
            elseif enum == this.ENUM_MODE_TXAMPFILE
                wts = this.amptx;
            else
                wts = ones(size(this.amp));
            end
            
                figure;
                hs=scatter3( ...
                    xData, ...
                    yData, ...
                    zData, ...
                    this.plot_marker_size,...
                    20*log10(abs(wts)),'filled');
                set(hs,'Marker','sq');
                colorbar;
                colormap jet;
               
            if     enum==this.ENUM_MODE_AMPFILE
                title('Rx Sum Weights Magnitude (dB)');
            elseif enum == this.ENUM_MODE_DAZ_AMPFILE
                title('Rx Dif AZ Weights Magnitude (dB)');
            elseif enum == this.ENUM_MODE_DEL_AMPFILE
                title('Rx Dif EL Weights Magnitude (dB)');
            elseif enum == this.ENUM_MODE_TXAMPFILE
                title('Tx Weights Magnitude (dB)');
            else
                title('Tx Weights (dB)');
            end
            
            if(this.upd_isvalid == 0)
                view(0,90);
                xlabel('X Direction (m)');
                ylabel('Y Direction (m)');
                zlabel('Z Direction (m)');
            else
                set(gca,'YDir','reverse','DataAspectRatio',[1 1 1]);
                ylabel('U (Y), meters');
                zlabel('V (Z), meters');
                xlabel('W (X), Normal, meters');
                set(hs,'MarkerEdgeColor',0.5*[1 1 1]);
                view(-90,0);
            end
            
            if ~isreal(wts)
                figure;
                hs=scatter3( ...
                    xData, ...
                    yData, ...
                    zData, ...
                    this.plot_marker_size,...
                    angle(wts)*180/pi,'filled');
                set(hs,'Marker','sq');
                colorbar;
                colormap jet;
                
                if     enum==this.ENUM_MODE_AMPFILE
                    title('Rx Sum Weights Phase (deg)');
                elseif enum == this.ENUM_MODE_DAZ_AMPFILE
                    title('Rx Dif AZ Weights Phase (deg)');
                elseif enum == this.ENUM_MODE_DEL_AMPFILE
                    title('Rx Dif EL Weights  Phase (deg)');
                elseif enum == this.ENUM_MODE_TXAMPFILE
                    title('Tx Weights  Phase (deg)');
                else
                    title('Tx Weights  Phase (deg)');
                end
                
                if(this.upd_isvalid == 0)
                    view(0,90);
                    xlabel('X Direction (m)');
                    ylabel('Y Direction (m)');
                    zlabel('Z Direction (m)');
                else
                    set(gca,'YDir','reverse','DataAspectRatio',[1 1 1]);
                    ylabel('U (Y), meters');
                    zlabel('V (Z), meters');
                    xlabel('W (X), Normal, meters');
                    set(hs,'MarkerEdgeColor',0.5*[1 1 1]);
                    view(-90,0);
                end
                
            end
            
            

            
        end
        
        function [g] = plotFullPattern(this,azvec,elvec,mode,space)
            if(~exist('azvec'))
                fprintf('Syntax: plotFullPattern(<ant>,[-90:1:90]*pi/180,[-90:1:90]*pi/180,<ant>.ENUM_MODE_???,space)\nspace can be either azel or sine');
                return;
            end
            if(~exist('space'))
                space = 'azel';
            end
            
            if(strcmpi(space,'azel'))
                g=getGainTable( ...
                    this,azvec,elvec, ...
                    mode);
            else
                %note: azvec and elvec here are u and v vectors
                g=getGainTableSinSpc( ...
                    this,azvec,elvec, ...
                    mode);
            end
            
            maxv = max(abs(g(:)));
            maxp = 20*log10(maxv);
            minp = maxp - 60;
            
            if(strcmpi(space,'azel'))
                hi=imagesc(180/pi*azvec,180/pi*elvec,20*log10(abs(g)));
                xlabel('Azimuth(deg)');
                ylabel('Elevation(deg)');
                grid on;
            else
                hi=imagesc(azvec,elvec,20*log10(abs(g)));
                xlabel('u');
                ylabel('v');
                set(hi,'AlphaData',(g~=0));
                grid on;
            end
            set(gca,'YDir','normal');
            
            if(strcmp(this.name,''))
                title('Antenna Pattern Power Response');
            else
                title({'Antenna Pattern Power Response',this.name});
            end
            colorbar;
            caxis([minp maxp]);
        end
        
        function [g] = upd_plotFullPattern(this,azvec,elvec,mode,space)
            if(~exist('azvec'))
                fprintf('Syntax: plotFullPattern(<ant>,[-90:1:90]*pi/180,[-90:1:90]*pi/180,<ant>.ENUM_MODE_???,space)\nspace can be either azel or sine');
                return;
            end
            if(~exist('space'))
                space = 'azel';
            end
            
            if(strcmpi(space,'azel'))
                g=upd_getGainTable( ...
                    this,azvec,elvec, ...
                    mode);
            else
                %note: azvec and elvec here are u and v vectors
                g=getGainTableSinSpc( ...
                    this,azvec,elvec, ...
                    mode);
            end
            
            maxv = max(abs(g(:)));
            maxp = 20*log10(maxv);
            minp = maxp - 60;
            
            if(strcmpi(space,'azel'))
                hi=imagesc(180/pi*azvec,180/pi*elvec,20*log10(abs(g)));
                xlabel('Azimuth(deg)');
                ylabel('Elevation(deg)');
                grid on;
            else
                hi=imagesc(azvec,elvec,20*log10(abs(g)));
                xlabel('u');
                ylabel('v');
                set(hi,'AlphaData',(g~=0));
                grid on;
            end
            set(gca,'YDir','normal');
            
            if(strcmp(this.name,''))
                title('Antenna Pattern Power Response');
            else
                title({'Antenna Pattern Power Response',this.name});
            end
            colorbar;
            caxis([minp maxp]);
        end
        
        function this = matlab_setElementPattern(this,angle_rad,gain_v)
            this.elemFacPresent = 1;
            if(angle_rad(1) ~= -pi/2 || angle_rad(end) ~= pi/2)
                error('setElementPattern: bounds must be +/- pi/2');
            end
            
            this.elementPattern = [angle_rad gain_v];
            
            % Populate updated array element patterns
            this.upd_elementPatternPresent = 1;
            [this.upd_elementPatternAzG,this.upd_elementPatternElG] = ...
                meshgrid(angle_rad,angle_rad);
            
            % Compute solid angle from 0,0
            solidAngle = acos(cos(this.upd_elementPatternAzG).*cos(this.upd_elementPatternElG));
            this.upd_elementPattern = interp1(angle_rad,gain_v,solidAngle);
            
            % Populate peak element gain as scalar value (for reference)
            this.upd_elementGaindBV = max(this.upd_elementPattern(:));
            
        end
        
        function this = matlab_setElementPattern2D(this,az_rad,el_rad,gain_v)
            this.elemFacPresent = 1;
%             this.elementPattern(:,1) = angle_rad;
%             this.elementPattern(:,2) = gain_v;
            
            % Populate updated array element patterns
            this.upd_elementPatternPresent = 1;
            [this.upd_elementPatternAzG,this.upd_elementPatternElG] = ...
                meshgrid(az_rad,el_rad);
            
            this.elemAz = this.upd_elementPatternAzG;
            this.elemEl = this.upd_elementPatternElG;
            this.elementPattern = gain_v;
            this.upd_elementPattern = gain_v;
            
            
            % Populate peak element gain as scalar value (for reference)
            this.upd_elementGaindBV = max(this.upd_elementPattern(:));
            
        end
        
        function gain = getGain(this,azInertial,elInertial,mode)
            if(this.upd_isvalid == 1)
                error('This array is using the updated storage & calculations. Use upd_getGain()');
            end
            
            [nRow,nCol] = size(azInertial);
            
            % if az/el inputs are vectors, ensure they are ROW VECTORS
            azInertial = reshape(azInertial,1,[]); % JAH - vectorized inputs
            elInertial = reshape(elInertial,1,[]); % JAH - vectorized inputs
            
            % Compute the antenna-relative angles based upon installation
            % Note: antenna is assumed to have its normal (z) aligned with
            % the host platform's longitudinal axis. This can be changed
            % using the changeInstNormal(az,el) method
            [az,el]=this.xformAzEl(azInertial,elInertial);
            
            if(~exist('mode','var'))
                mode = this.ENUM_MODE_UNIFORM;
            end
            
            % convert az,el (in antenna frame) to unit vector
            [uz,ux,uy]=sph2cart(az,el,1);
            %u = [ux uy uz]; % JAH - original line
            u = [ux.' uy.' uz.']; % JAH - vectorized inputs
            
            % NOTE: z is normal
            % Compute total steering angle
            incAngle = acos(uz);
            
            if(this.needPhaseQuantized == 1)
                applyPhaseQuantization(this)
            end
            
            % compute element factor for this far-field location
            if(this.elemFacPresent == 1)
                %assuming here that the importFEKOelement routine was used
                %need to set this up to only run once - CTM, clayton
                if(any([nRow nCol] == 1))
                    [azMesh,elMesh] = meshgrid(az,el);
                else
                    azMesh = reshape(az,nRow,nCol);
                    elMesh = reshape(el,nRow,nCol);
                end
                
                try
                    this_elFac_real = interp2(this.elemAz,this.elemEl,real(this.elementPattern),azMesh,elMesh,'nearest');
                    this_elFac_imag = interp2(this.elemAz,this.elemEl,imag(this.elementPattern),azMesh,elMesh,'nearest');
                    this_elFac = this_elFac_real + 1i*this_elFac_imag;
                    this_elFac = reshape(this_elFac,1,[]);
                catch
                    warning('importFEKOelement assumption broken. using 1-D lookup');
                    this_elFac = interp1(this.elementPattern(:,1),this.elementPattern(:,2),incAngle);
                end
            else
                % apply standard cos(theta) scanning loss for isotropic radiator assumption
                this_elFac = sqrt(abs(cos(incAngle).^this.cos_pow_scn));
            end
            csum = 0;
            
            if(this.calcMethod == 1) % implemented as fortran does
                
                for k = 1:this.nElements
                    this_a = 1 + this.aVar * randn();
                    this_p = 0 + this.pVar * randn();
                    
                    this.phs(k) = (ux*this.pos(1,k) + ...
                        uy*this.pos(2,k) + ...
                        uz*this.pos(3,k)) * this.twopi_ovr_lambda + ...
                        this_p;
                    this.phs(k) = this.phs(k) - this.phs_steer(k);
                    
                    val = exp(1i*this.phs(k));
                    
                    if(mode == this.ENUM_MODE_AMPFILE)
                        csum = csum + (this_a * this.amp(k)) * val;
                    elseif(mode == this.ENUM_MODE_UNIFORM)
                        csum = csum + (this_a              ) * val;
                    elseif(mode == this.ENUM_MODE_DAZ_AMPFILE)
                        csum = csum + (this_a * this.ampaz(k)) * val;
                    elseif(mode == this.ENUM_MODE_DEL_AMPFILE)
                        csum = csum + (this_a * this.ampel(k)) * val;
                    end
                end
                
            else % vectorized for matlab
                % summation is over columns
                this_a = 1 + this.aVar * randn(1,this.nElements,this.precisionClass);
                this_p = 0 + this.pVar * randn(1,this.nElements,this.precisionClass);
                %                 this_a = this_a .* (~this.element_off);
                %this.phs = ( ...
                %    ux * this.pos(1,:) + ...
                %    uy * this.pos(2,:) + ...
                %    uz * this.pos(3,:)) * this.twopi_ovr_lambda + this_p - this.phs_steer;
                
                %this.phs = (u * this.pos) * this.twopi_ovr_lambda + ...
                %    this_p - this.phs_steer + this.phs_sphcmp;  % JAH - original
                
                onecolvec = ones(size(u,1),1);
                
                this.phs = (u * this.pos) * this.twopi_ovr_lambda + ...
                    onecolvec*(this_p - this.phs_steer + this.phs_sphcmp); % JAH - vectorized inputs
                
                val = exp(1i*this.phs);
                
                if(mode == this.ENUM_MODE_AMPFILE)
                    csum = onecolvec*(this_a .* this.amp) .* val;
                elseif(mode == this.ENUM_MODE_UNIFORM)
                    csum = onecolvec*(this_a              ) .* val;
                elseif(mode == this.ENUM_MODE_DAZ_AMPFILE)
                    csum = onecolvec*(this_a .* this.ampaz) .* val;
                elseif(mode == this.ENUM_MODE_DEL_AMPFILE)
                    csum = onecolvec*(this_a .* this.ampel) .* val;
                end
                %csum = sum(csum); % JAH - original
                csum = sum(csum,2).*this_elFac.'; % JAH - vectorized inputs, 2nd arg sum over N-elements dimensions (sum each row)
                csum = reshape(csum,1,[]); % JAH - vectorized inputs, put back into ROW VECTOR
            end
            
            % jcl(-)gain = csum * this.neScale * this.efficiency * this.vPeak * this_elFac;
            
            % KDS - this if block no longer needed due to proper
            % normalization below (/sqrt(#elements))
            %if(this.elemFacPresent == 1)
            %    gain = csum .* this_elFac; % jcl(+)
            %else
            %    gain = csum .* this.neScale .* this.vPeak .* this_elFac; % jcl(+)
            %end
            
            gain = this_elFac *  csum ./ sqrt(this.nElements);
            
            %             AF = sqrt(csum);
            %             EF = this_elFac;
            %             gain = AF .* EF .* this.efficiency;
            
            %Put back into same size as input CTM
            gain = reshape(gain,nRow,nCol);
        end
        
        function setSteer(this,theta,phi,idx_elements)
            if(this.upd_isvalid == 1)
                error('This array is using the updated storage & calculations. Use upd_setSteer()');
            end
            
            this.azSteer = theta;
            this.elSteer = phi;
            [uz,ux,uy]=sph2cart(theta,phi,1.0);
            if(~exist('idx_elements'))
                idx_elements = 1 : this.nElements;
            end
            for k = idx_elements
                this.phs_steer(k) =( ...
                    ux*this.pos(1,k) + ...
                    uy*this.pos(2,k) + ...
                    uz*this.pos(3,k) ) * this.twopi_ovr_lambda;
                %this.phs_steer(k) = this%phs_steer(k)
            end
            if(this.pShiftQuantization == 1)
                this.needPhaseQuantized = 1;
                applyPhaseQuantization(this)
            end
        end
        
        %If this is called, it will set flags that will cause phase steering to be
        %quantized via applyPhaseQuantization from setSteer.
        function setPhaseQuantization(this,minPhaseDeg,maxPhaseDeg,nBits)
            if(this.upd_isvalid == 1)
                error('Updated array. Use upd_setPhaseQuantBits');
            end
            this.pBits              = nBits;
            this.pMin               = minPhaseDeg * pi / 180;
            this.pMax               = maxPhaseDeg * pi / 180;
            this.pStep              = (this.pMax - this.pMin) ...
                / 2^(this.pBits-1);
            this.pShiftQuantization = 1;
            this.needPhaseQuantized = 1;
        end
        
        function upd_setPhaseQuantBits(this,nBits)
            if(nBits > 0)
                this.pBits              = nBits;
                this.pStep              = 2*pi ./ 2^this.pBits;
                this.pShiftQuantization = 1;
                this.needPhaseQuantized = 1;
            else
                this.pShiftQuantization = 0;
                this.needPhaseQuantized = 0;
            end
        end
        
        %quantize the phase according to values set by
        %setPhaseQuantization
        function applyPhaseQuantization(this)
            if(this.upd_isvalid == 1)
                error('Updated array. Use upd_applyPhaseQuantization');
            end
            for k = 1:this.nElements
                this.phs_steer(k) = this.pStep * floor( ...
                    (this.phs_steer(k) - this.pMin + .5*this.pStep) ...
                    / this.pStep );
            end
            this.needPhaseQuantized = 0;
        end
        
        function upd_applyPhaseQuantization(this)
            this.phs_steer = this.pStep * floor(mod(this.phs_steer,2*pi) ./ (this.pStep));
            this.needPhaseQuantized = 0;
        end
        
        function [scaleFac_az,scaleFac_el]=angDiscCal(this,az,el,aerng,aestep)
            nSteps = aerng/aestep + 1;
            
            % 1 - steer beam to an az el
            this.setSteer(az,el);
            
            beta_N = 0.0;
            beta_D = 0.0;
            
            
            % For all angle steps around the discriminator cal
            for k = 1:nSteps
                % compute this angle lookup
                azi = az - 0.5 * aerng + (k-1)*aestep; % should span from steered az -/+ aerng
                sum_v = getGain(this,azi,el,this.ENUM_MODE_AMPFILE);
                del_v = getGain(this,azi,el,  this.ENUM_MODE_DAZ_AMPFILE);
                
                sum_v = abs(sum_v);
                del_vi = imag(del_v);
                del_vr = real(del_v);
                %                 hold on;
                %                 plot(azi,[sum_v del_vi del_vr],'.');
                
                
                discr = del_vi / sum_v;
                
                beta_N = beta_N + discr * (azi - az);
                beta_D = beta_D + (azi-az) * (azi-az);
            end
            scaleFac_az = beta_N / beta_D;
            
            beta_N = 0.0;
            beta_D = 0.0;
            
            % For all angle steps around the discriminator cal
            for k = 1:nSteps
                % compute this angle lookup
                eli = el - 0.5 * aerng + (k-1)*aestep; % should span from steered az -/+ aerng
                sum_v = getGain(this,az,eli,  this.ENUM_MODE_AMPFILE);
                del_v = getGain(this,az,eli,  this.ENUM_MODE_DEL_AMPFILE);
                %write(701) sum_v,Del_v,az,eli
                sum_v = abs(sum_v);
                del_v = imag(del_v);
                
                discr = del_v / sum_v;
                
                %write(700) 2.0d0,el,eli,sum_v,Del_v,discr
                
                beta_N = beta_N + discr * (eli - el);
                beta_D = beta_D + (eli - el) * (eli - el);
            end
            scaleFac_el = beta_N / beta_D;
        end
        
        function res = getGainTable(this,azpts,elpts,mode)
            tic;
            res = zeros(length(elpts),length(azpts),this.precisionClass);
            hw = waitbar(0,sprintf('Computing Far-Field Gain: %d%% Complete',0));
            for kaz = 1:length(azpts)
                for kel = 1:length(elpts)
                    res(kel,kaz) = getGain(this,azpts(kaz),elpts(kel),mode);
                end
                waitbar(kaz/length(azpts),hw,sprintf('Computing Far-Field Gain: %d%% Complete',floor(100*kaz/length(azpts))));
            end
            close(hw);
            tout = toc;
            fprintf('Total Time: %f sec\nTime per Az/El: %f sec\n',tout,tout/(length(azpts)*length(elpts)));
        end
        
        function res = upd_getGainTable(this,azpts,elpts,mode)
            tic;
            res = zeros(length(elpts),length(azpts),this.precisionClass);
            hw = waitbar(0,sprintf('Computing Far-Field Gain: %d%% Complete',0));
            for kaz = 1:length(azpts)
                for kel = 1:length(elpts)
                    res(kel,kaz) = upd_getGain(this,azpts(kaz),elpts(kel),mode);
                end
                waitbar(kaz/length(azpts),hw,sprintf('Computing Far-Field Gain: %d%% Complete',floor(100*kaz/length(azpts))));
            end
            close(hw);
            tout = toc;
            fprintf('Total Time: %f sec\nTime per Az/El: %f sec\n',tout,tout/(length(azpts)*length(elpts)));
        end
        
        function res = getSupGainTable(this,azpts,elpts,mode)
            res = zeros(length(elpts),length(azpts),this.precisionClass);
            for kaz = 1:length(azpts)
                for kel = 1:length(elpts)
                    res(kel,kaz) = getGain(this,azpts(kaz),elpts(kel),mode);
                end
            end
        end
        
        function res = getGainTableSinSpc(this,upts,vpts,mode)
            tic;
            res = zeros(length(vpts),length(upts),this.precisionClass);
            hw = waitbar(0,sprintf('Computing Far-Field Gain: %d%% Complete',0));
            for ku = 1:length(upts)
                for kv = 1:length(vpts)
                    w = sqrt(1 - upts(ku)^2 - vpts(kv)^2);
                    ir = isreal(w);
                    w = real(w);
                    
                    %az = asin(w);
                    %el = atan2(vpts(kv),upts(ku));
                    if(ir)
                        [az,el,r]=cart2sph(w,-vpts(kv),-upts(ku));
                        res(ku,kv) = getGain(this,az,el,mode);
                    else
                        res(ku,kv)=0;
                    end
                    %res(ku,kv) = complex(az,el);
                end
                waitbar(ku/length(upts),hw,sprintf('Computing Far-Field Gain: %d%% Complete',floor(100*ku/length(upts))));
            end
            close(hw);
            tout = toc;
            fprintf('Total Time: %f sec\nTime per Az/El: %f sec\n',tout,tout/(length(upts)*length(vpts)));
        end
        
        function [resaz,resel] = getDiscrTables(this,azpts,elpts,angbound,angstep)
            % all are in radians
            resaz = zeros(length(elpts),length(azpts),this.precisionClass);
            resel = zeros(length(elpts),length(azpts),this.precisionClass);
            hw = waitbar(0,sprintf('Computing Discriminator Slopes: %d%% Complete',0));
            for kaz = 1:length(azpts)
                for kel = 1:length(elpts)
                    [af,ef]=this.angDiscCal(azpts(kaz),elpts(kel), ...
                        angbound,angstep);
                    resaz(kel,kaz) = af;
                    resel(kel,kaz) = ef;
                end
                waitbar(kaz/length(azpts),hw,sprintf('Computing Discriminator Slopes: %d%% Complete',floor(100*kaz/length(azpts))));
            end
            close(hw);
            tout = toc;
            fprintf('Total Time: %f sec\nTime per Az/El: %f sec\n',tout,tout/(length(azpts)*length(elpts)));
        end
        
        function res = getAzCut(this,azVec,atElevation,mode,outputStyle)
            % getAzCut(this,azVec,atElevation,mode,outputStyle)
            % all in radians, output style either 'P' power or 'CV' complex
            % voltage
            if(~exist('outputStyle'))
                outputStyle = 'CV';
            end
            res = zeros(1,length(azVec),this.precisionClass);
            for kaz = 1:length(azVec)
                res(kaz) = getGain(this,azVec(kaz),atElevation,mode);
            end
            if(strcmp(outputStyle,'P'))
                res = 20*log10(abs(res));
            end
        end
        
        function res = upd_getAzCut(this,azVec,atElevation,mode,outputStyle)
            % getAzCut(this,azVec,atElevation,mode,outputStyle)
            % all in radians, output style either 'P' power or 'CV' complex
            % voltage
            if(~exist('outputStyle'))
                outputStyle = 'CV';
            end
            res = zeros(1,length(azVec),this.precisionClass);
            for kaz = 1:length(azVec)
                res(kaz) = upd_getGain(this,azVec(kaz),atElevation,mode);
            end
            if(strcmp(outputStyle,'P'))
                res = 20*log10(abs(res));
            end
        end
        
        function res = getElCut(this,elVec,atAzimuth,mode,outputStyle)
            if(~exist('outputStyle'))
                outputStyle = 'CV';
            end
            res = zeros(1,length(elVec),this.precisionClass);
            for kel = 1:length(elVec)
                res(kel) = getGain(this,atAzimuth,elVec(kel),mode);
            end
            if(strcmp(outputStyle,'P'))
                res = 20*log10(abs(res));
            end
        end
        
        function res = upd_getElCut(this,elVec,atAzimuth,mode,outputStyle)
            if(~exist('outputStyle'))
                outputStyle = 'CV';
            end
            res = zeros(1,length(elVec),this.precisionClass);
            for kel = 1:length(elVec)
                res(kel) = upd_getGain(this,atAzimuth,elVec(kel),mode);
            end
            if(strcmp(outputStyle,'P'))
                res = 20*log10(abs(res));
            end
        end
        
        
        
        function compTaylorWgt(this,NBAR,SLL)
            % Syntax: compTaylorWgt(this,NBAR,SLL)
            % this requires access to MATLAB taylorwin function
            % SLL is input as a negative [dB]
            % compute a linear taylor of length = # elements
            taylor_ref_x = linspace(-1,1,this.nElements);
            taylor_ref_y = taylorwin(this.nElements,NBAR,SLL);
            
            % normalize taylor output
            %taylor_ref_y = taylor_ref_y ./ max(abs(taylor_ref_y(:)));
            
            
            % compute normalized range (0 = middle, +/- 1 = edge)
            if(~isempty(this.rng))
                nrng = this.rng ./ max(abs(this.rng(:)));
            else
                errordlg('rng property not defined! Can''t build taylor weighting!','Error');
                error('rng property not defined!');
            end
            
            % lookup amplitude weighting based upon range
            this.amp = interp1(taylor_ref_x,taylor_ref_y,nrng);
        end
        
        function compUniformWgt(this)
            this.amp = 0*this.pos(1,:) + 1;
        end
        
        function compBaylissWgt(this,NBAR,SLL)
            % compBaylissWgt(this,NBAR,SLL)
            %
            wba = GTRI_cbayliss([this.pos(1,:)' this.pos(2,:)'],SLL,NBAR).';
            wbe = GTRI_cbayliss([this.pos(2,:)' this.pos(1,:)'],SLL,NBAR).';
            wba = sign(wba) .* sqrt(abs(wba));
            wbe = sign(wbe) .* sqrt(abs(wbe));
            %wba(isnan(re)) = nan;
            %wbe(isnan(re)) = nan;
            this.ampaz = wba;
            this.ampel = wbe;
        end
        
        function cumPower = computeCumulativeNoise(this,mode)
            % Purpose: Compute cumulative noise power when digitally
            % scaling analog noise on each element to achieve taper
            
            nAvg = 1000;
            cumPower = 0;
            for k = 1 : nAvg
                if(this.upd_isvalid == 1)
                    nzSig = randn(this.nElements,1) + 1i.*randn(this.nElements,1);
                else
                    nzSig = randn(1,this.nElements) + 1i.*randn(1,this.nElements);
                end
                
                switch(mode)
                    case this.ENUM_MODE_UNIFORM
                        nzSig = nzSig .* 1.0;
                        
                    case this.ENUM_MODE_AMPFILE
                        nzSig = nzSig .* this.amp;
                        
                    case this.ENUM_MODE_DAZ_AMPFILE
                        nzSig = nzSig .* this.ampaz;
                        
                    case this.ENUM_MODE_DEL_AMPFILE
                        nzSig = nzSig .* this.ampel;
                end
                
                cumPower = cumPower + abs(sqrt(mean(abs(nzSig).^2))).^2;
            end
            cumPower = cumPower ./ nAvg;
            
        end
        
        function writeFortranFiles(this,arrayName)
            fid = fopen([arrayName '_Elements.txt'],'w');
            fprintf(fid,'%d\n',this.nElements);
            for k = 1 : this.nElements
                fprintf(fid,'%18.6f %18.6f %18.6f\n',this.pos(:,k));
            end
            fclose(fid);
            
            fid = fopen([arrayName '_AmpWeight.txt'],'w');
            fprintf(fid,'%d\n',this.nElements);
            for k = 1 : this.nElements
                fprintf(fid,'%18.6f\n',this.amp(k));
            end
            fclose(fid);
            
            fid = fopen([arrayName '_AmpAzWeight.txt'],'w');
            fprintf(fid,'%d\n',this.nElements);
            for k = 1 : this.nElements
                fprintf(fid,'%18.6f\n',this.ampaz(k));
            end
            fclose(fid);
            
            fid = fopen([arrayName '_AmpElWeight.txt'],'w');
            fprintf(fid,'%d\n',this.nElements);
            for k = 1 : this.nElements
                fprintf(fid,'%18.6f\n',this.ampel(k));
            end
            fclose(fid);
        end
        
        function plotSuite(this)
            figure;
            this.plotGeometry();
            
            figure;
            this.plotGeometry(this.amp);
            title('Weighted Amplitude (mag)');
            
            figure;
            this.plotGeometry(180/pi*this.phs_steer);
            title('Phasing (steering,deg)');
            
            figure;
            this.plotGeometry(this.ampaz);
            title('Weighted Amplitude (\DeltaAz)');
            
            figure;
            this.plotGeometry(this.ampel);
            title('Weighted Amplitude (\DeltaEl)');
            
            % Cuts
            azd = [-90:this.plot_2d_res_deg:90];
            eld = [-90:this.plot_2d_res_deg:90];
            azdr = pi/180*[-90:this.plot_2d_res_deg:90];
            eldr = pi/180*[-90:this.plot_2d_res_deg:90];
            azcu = this.getAzCut(azdr,0.0,this.ENUM_MODE_UNIFORM,'P');
            azcw = this.getAzCut(azdr,0.0,this.ENUM_MODE_AMPFILE,'P');
            azcd = this.getAzCut(azdr,0.0,this.ENUM_MODE_DAZ_AMPFILE,'P');
            elcu = this.getElCut(eldr,0.0,this.ENUM_MODE_UNIFORM,'P');
            elcw = this.getElCut(eldr,0.0,this.ENUM_MODE_AMPFILE,'P');
            elcd = this.getElCut(eldr,0.0,this.ENUM_MODE_DEL_AMPFILE,'P');
            
            figure;
            plot(azd,[azcu.' azcw.' azcd.']);
            xlabel('Azimuth (deg)');
            ylabel('Gain (dBi)');
            grid on;
            title('Azimuth Cut');
            legend('Uniform','Weighted','Diff');
            
            figure;
            plot(eld,[elcu.' elcw.' elcd.']);
            xlabel('Elevation (deg)');
            ylabel('Gain (dBi)');
            grid on;
            title('Elevation Cut');
            legend('Uniform','Weighted','Diff');
            
            figure;
            plotFullPattern(this,azdr,eldr,this.ENUM_MODE_UNIFORM);
            title('Power Gain (dBi), Uniform');
            
            figure;
            plotFullPattern(this,azdr,eldr,this.ENUM_MODE_AMPFILE);
            title('Power Gain (dBi), Weighted');
            
            % jcl(-)            figure;
            % jcl(-)            plotFullPattern(this,azdr,eldr,this.ENUM_MODE_DAZ_AMPFILE);
            % jcl(-)            title('Power Gain (dBi), \Delta Az');
            
            % jcl(-)            figure;
            % jcl(-)            plotFullPattern(this,azdr,eldr,this.ENUM_MODE_DEL_AMPFILE);
            % jcl(-)            title('Power Gain (dBi), \Delta El');
            
        end
        
        function upd_plotSuite(this)
            figure;
            this.plotGeometry();
            
            figure;
            this.plotGeometry(this.amp);
            title('Weighted Amplitude (mag)');
            
            figure;
            this.plotGeometry(180/pi*this.phs_steer);
            title('Phasing (steering,deg)');
            
            figure;
            this.plotGeometry(this.ampaz);
            title('Weighted Amplitude (\DeltaAz)');
            
            figure;
            this.plotGeometry(this.ampel);
            title('Weighted Amplitude (\DeltaEl)');
            
            % Cuts
            azd = [-90:this.plot_2d_res_deg:90];
            eld = [-90:this.plot_2d_res_deg:90];
            azdr = pi/180*[-90:this.plot_2d_res_deg:90];
            eldr = pi/180*[-90:this.plot_2d_res_deg:90];
            azcu = this.upd_getAzCut(azdr,0.0,this.ENUM_MODE_UNIFORM,'P');
            azcw = this.upd_getAzCut(azdr,0.0,this.ENUM_MODE_AMPFILE,'P');
            azcd = this.upd_getAzCut(azdr,0.0,this.ENUM_MODE_DAZ_AMPFILE,'P');
            elcu = this.upd_getElCut(eldr,0.0,this.ENUM_MODE_UNIFORM,'P');
            elcw = this.upd_getElCut(eldr,0.0,this.ENUM_MODE_AMPFILE,'P');
            elcd = this.upd_getElCut(eldr,0.0,this.ENUM_MODE_DEL_AMPFILE,'P');
            
            figure;
            plot(azd,[azcu.' azcw.' azcd.']);
            xlabel('Azimuth (deg)');
            ylabel('Gain (dBi)');
            grid on;
            title('Azimuth Cut');
            legend('Uniform','Weighted','Diff');
            
            figure;
            plot(eld,[elcu.' elcw.' elcd.']);
            xlabel('Elevation (deg)');
            ylabel('Gain (dBi)');
            grid on;
            title('Elevation Cut');
            legend('Uniform','Weighted','Diff');
            
            figure;
            upd_plotFullPattern(this,azdr,eldr,this.ENUM_MODE_UNIFORM);
            title('Power Gain (dBi), Uniform');
            
            figure;
            upd_plotFullPattern(this,azdr,eldr,this.ENUM_MODE_AMPFILE);
            title('Power Gain (dBi), Weighted');
            
            % jcl(-)            figure;
            % jcl(-)            plotFullPattern(this,azdr,eldr,this.ENUM_MODE_DAZ_AMPFILE);
            % jcl(-)            title('Power Gain (dBi), \Delta Az');
            
            % jcl(-)            figure;
            % jcl(-)            plotFullPattern(this,azdr,eldr,this.ENUM_MODE_DEL_AMPFILE);
            % jcl(-)            title('Power Gain (dBi), \Delta El');
            
        end
        
        
        % jcl(+)->
        function bw = get3dbBeamWidth(this,cut_angles,cut_gain)
            % cut_angles vector in degrees
            % cut_gain vector in dBi
            max_cut_gain_less3db = max(cut_gain) - 3;
            bw =( max(find(cut_gain >= max_cut_gain_less3db)) - min(find(cut_gain >= max_cut_gain_less3db)) ) * (cut_angles(end) - cut_angles(max(end - 1,1)));
        end
        
        function writeSupAntPat(this,...        % Array data structure
                DataLbl,...                     % Array description, string
                OutputDataPath,...              % Directory for saved pattern data, string
                SupAntFilNam,...                % Suppressor ANTENNA-PATTERN file name with file extension (e.g. ant_pat.dat, or Z://mypatterns/array.ant etc.), string
                ClassLbl,...                    % Classification desription of data, string
                FreqIntrvl_GHz,...              % Antenna pattern frequency intervals (ascending order), real array
                FreqVal_GHz,...                 % Frequency values corresponding to each frequency interval (ascending order), real array
                AzIntrvl_deg,...                % Azimuth intervals (ascending order, i.e. from right to left), real array
                AzVal_deg,...                   % Azimuth values corresponding to each azimuth interval (ascending order), real array
                AzSteerAngle_deg,...            % Azimuth mainbeam steering angle (deg), real value
                ElIntrvl_deg,...                % Elevation intervals (ascending order, i.e. from bottom to top), real array
                ElVal_deg,...                   % Elevation values corresponding to each elevation interval (ascending order), real array
                ElSteerAngle_deg,...            % Elevation mainbeam steering angle (deg), real value
                ElementXYZpos_m,...             % Element x,y,z positions (m), real array
                HpolElementArea_sqm,...         % Element area (square-meters) for horizontal polarization corresponding to each frequency value, real array
                VpolElementArea_sqm,...         % Element area (square-meters) for vertical polarization corresponding to each frequency value, real array
                ElementAmp,...                  % Element amplitudes (weights), real array
                ArrayEff,...                    % Overall array efficiency [0~1], real value
                AmpVar_V,...                    % Amplitude variation over array (+/- Volts)
                PhaseVar_deg,...                % Phase variation over array (+/- deg)
                HpolElementGainsForFreq_dBi,... % Horizontal polarization elemental gains (dBi) corresponding to each frequency value, real array
                VpolElementGainsForFreq_dBi...  % Vertical polarization elemental gains (dBi) corresponding to each frequency value, real array
                )
            
            %----------------------------------------------------------------------------------
            % Note, the output statements below are formatted with the following assumptions:
            %   1.  HH-Pol & VV-Pol (Tx/Rx) polarizations are always output to the file.
            %       If only one polarization is available, make them duplicates.
            %   2.  Output data units are:
            %         GHZ(frequency)
            %         DEG(angle)
            %         DB(gain)
            %   3.  Frequency output precision is 1 MHz, frequency range=[1~99999] MHz
            %   4.  Azimuth output precision is 0.1 deg, azimuth range=[-180~180] deg
            %   5.  Elevation output precision is 0.1 deg, azimuth range=[-90~90] deg
            %   6.  Gain output precision is 0.01 dB, gain range=[-999.99~999.99] dB
            %   7.  Suppressor files are limited to 72 columns
            %----------------------------------------------------------------------------------
            % Initialize gain data structure
            Gain_dB = zeros(length(ElVal_deg),length(AzVal_deg));
            % Open a file with the user specified name
            FullSupAntFilNam = strcat(OutputDataPath,SupAntFilNam,'.dat');
            fid = fopen([FullSupAntFilNam],'w');
            % Output Suppressor ANTENNA-PATTERN table header
            fprintf(fid,'ANTENNA-PATTERN\n');
            % Ouput classification label
            fprintf(fid,'$ %s\n',ClassLbl);
            % Ouput data label
            fprintf(fid,'$ %s\n',DataLbl);
            fprintf(fid,' DIMENSION 1  POL\n');
            fprintf(fid,'  HORIZ-POL  VERT-POL\n');
            % Output HH-Pol (Tx/Rx-Polarization) data
            fprintf(fid,'  $ Freq Data for HHpol (Tx/Rx-Polarization)\n');
            fprintf(fid,'  DIMENSION 2  FREQ (GHZ)\n');
            fprintf(fid,'   %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',FreqIntrvl_GHz); % Frequency output precision is 1 MHz
            fprintf(fid,'\n'); % Ensure start of a new line
            for nF = 1:length(FreqIntrvl_GHz)-1
                fprintf(fid,'   $ Az Data for HHpol, F(GHz)=%6.3f[%6.3f~%6.3f]\n',FreqVal_GHz(nF),FreqIntrvl_GHz(nF),FreqIntrvl_GHz(nF+1));
                fprintf(fid,'   DIMENSION 3  AZ (DEG)\n');
                fprintf(fid,'    %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n',AzIntrvl_deg); % Azimuth output precision is 0.1 deg
                fprintf(fid,'\n'); % Ensure start of a new line
                % Array Configuration
                Fop_Hz = FreqVal_GHz(nF)*1e9;                                                                              % Array operating frequency (Hz)
                arr = array();                                                                                             % Create a new array object
                arr.setElementXYZ(ElementXYZpos_m(:,1),ElementXYZpos_m(:,2),ElementXYZpos_m(:,3),HpolElementArea_sqm(nF)); % Insert data into the array object
                arr.amp = ElementAmp;                                                                                      % Transpose input data from file
                Gelement_dBi = HpolElementGainsForFreq_dBi(nF);                                                            % Hpol element gain corresponding to input frequency
                ScaledArrayEff = 10^(Gelement_dBi/10)*ArrayEff;                                                            % Scaled array efficiency (applies measured element gain)
                arr.setTxProperties(Fop_Hz, ScaledArrayEff);                                                               % Efficiency (@ Operating Frequency)
                arr.aVar = AmpVar_V;                                                                                       % Amplitude variation over array (+/- V)
                arr.pVar = PhaseVar_deg*pi/180;                                                                            % Phase variation over array (+/- rad)
                arr.turnOn('all');                                                                                         % Turn on desired elements: 'all', 'none', or element #
                arr.setSteer(AzSteerAngle_deg/180*pi, ElSteerAngle_deg/180*pi);                                            % Az & El steering angles (rad)
                % Output progress to Matlab command window
                fprintf('Computing HH-pol gain pattern at %f GHz\n',FreqVal_GHz(nF));
                % Get gain table data
                gain2d = arr.getSupGainTable(AzVal_deg*pi/180,ElVal_deg*pi/180,arr.amp);
                % Convert gain from absolute voltage to power dBI
                Gain_dB = 20*log10(abs(gain2d));
                for nAz = 1:length(AzIntrvl_deg)-1
                    fprintf(fid,'    $ El Data for HHpol, F(GHz)=%6.3f[%6.3f~%6.3f], Az(deg)=%6.1f[%6.1f~%6.1f]\n',FreqVal_GHz(nF),FreqIntrvl_GHz(nF),FreqIntrvl_GHz(nF+1),AzVal_deg(nAz),AzIntrvl_deg(nAz),AzIntrvl_deg(nAz+1));
                    fprintf(fid,'    DIMENSION 4  EL (DEG)\n');
                    fprintf(fid,'     %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',ElIntrvl_deg); % Elevation output precision is 0.1 deg
                    fprintf(fid,'\n'); % Ensure start of a new line
                    fprintf(fid,'     $ Gain Data for HHpol, F(GHz)=%6.3f[%6.3f~%6.3f], Az(deg)=%6.1f[%6.1f~%6.1f], EL(deg)=[%5.1f~%5.1f]\n',FreqVal_GHz(nF),FreqIntrvl_GHz(nF),FreqIntrvl_GHz(nF+1),AzVal_deg(nAz),AzIntrvl_deg(nAz),AzIntrvl_deg(nAz+1),ElIntrvl_deg(1),ElIntrvl_deg(end));
                    fprintf(fid,'     GAIN (DB)\n');
                    fprintf(fid,'      %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',Gain_dB(:,nAz)); % Elevation output precision is 0.1 deg
                    fprintf(fid,'\n'); % Ensure start of a new line
                end;
            end;
            % Output VV-Pol (Tx/Rx-Polarization) data
            fprintf(fid,'  $ Freq Data for VVpol (Tx/Rx-Polarization)\n');
            fprintf(fid,'  DIMENSION 2  FREQ (GHZ)\n');
            fprintf(fid,'   %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',FreqIntrvl_GHz); % Frequency output precision is 1 MHz
            fprintf(fid,'\n'); % Ensure start of a new line
            for nF = 1:length(FreqIntrvl_GHz)-1
                fprintf(fid,'   $ Az Data for VVpol, F(GHz)=%6.3f[%6.3f~%6.3f]\n',FreqVal_GHz(nF),FreqIntrvl_GHz(nF),FreqIntrvl_GHz(nF+1));
                fprintf(fid,'   DIMENSION 3  AZ (DEG)\n');
                fprintf(fid,'    %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n',AzIntrvl_deg); % Azimuth output precision is 0.1 deg
                fprintf(fid,'\n'); % Ensure start of a new line
                % Array Configuration
                Fop_Hz = FreqVal_GHz(nF)*1e9;                                                                              % Array operating frequency (Hz)
                arr = array();                                                                                             % Create a new array object
                arr.setElementXYZ(ElementXYZpos_m(:,1),ElementXYZpos_m(:,2),ElementXYZpos_m(:,3),VpolElementArea_sqm(nF)); % Insert data into the array object
                arr.amp = ElementAmp;                                                                                      % Transpose input data from file
                Gelement_dBi = VpolElementGainsForFreq_dBi(nF);                                                            % Vpol element gain corresponding to input frequency
                ScaledArrayEff = 10^(Gelement_dBi/10)*ArrayEff;                                                            % Scaled array efficiency (applies measured element gain)
                arr.setTxProperties(Fop_Hz, ScaledArrayEff);                                                               % Efficiency (@ Operating Frequency)
                arr.aVar = AmpVar_V;                                                                                       % Amplitude variation over array (+/- V)
                arr.pVar = PhaseVar_deg*pi/180;                                                                            % Phase variation over array (+/- rad)
                arr.turnOn('all');                                                                                         % Turn on desired elements: 'all', 'none', or element #
                arr.setSteer(AzSteerAngle_deg/180*pi, ElSteerAngle_deg/180*pi);                                            % Az & El steering angles (rad)
                % Output progress to Matlab command window
                fprintf('Computing VV-pol gain pattern at %f GHz\n',FreqVal_GHz(nF));
                % Get gain table data
                gain2d = arr.getSupGainTable(AzVal_deg*pi/180,ElVal_deg*pi/180,arr.amp);
                % Convert gain from absolute voltage to power dBI
                Gain_dB = 20*log10(abs(gain2d));
                for nAz = 1:length(AzIntrvl_deg)-1
                    fprintf(fid,'    $ El Data for VVpol, F(GHz)=%6.3f[%6.3f~%6.3f], Az(deg)=%6.1f[%6.1f~%6.1f]\n',FreqVal_GHz(nF),FreqIntrvl_GHz(nF),FreqIntrvl_GHz(nF+1),AzVal_deg(nAz),AzIntrvl_deg(nAz),AzIntrvl_deg(nAz+1));
                    fprintf(fid,'    DIMENSION 4  EL (DEG)\n');
                    fprintf(fid,'     %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',ElIntrvl_deg); % Elevation output precision is 0.1 deg
                    fprintf(fid,'\n'); % Ensure start of a new line
                    fprintf(fid,'     $ Gain Data for VVpol, F(GHz)=%6.3f[%6.3f~%6.3f], Az(deg)=%6.1f[%6.1f~%6.1f], EL(deg)=[%5.1f~%5.1f]\n',FreqVal_GHz(nF),FreqIntrvl_GHz(nF),FreqIntrvl_GHz(nF+1),AzVal_deg(nAz),AzIntrvl_deg(nAz),AzIntrvl_deg(nAz+1),ElIntrvl_deg(1),ElIntrvl_deg(end));
                    fprintf(fid,'     GAIN (DB)\n');
                    fprintf(fid,'      %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',Gain_dB(:,nAz)); % Elevation output precision is 0.1 deg
                    fprintf(fid,'\n'); % Ensure start of a new line
                end;
            end;
            % Ouput classification label
            fprintf(fid,'$ %s\n',ClassLbl);
            % Output Suppressor ANTENNA-PATTERN table footer
            fprintf(fid,'END ANTENNA-PATTERN\n');
            % Close the user specified file
            fclose(fid);
        end
        % <-jcl(+)
        
        
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % This section contains KDS updates to the array code based on
        % 7/2019 deep dive for adaptive beamforming and MUSIC algorithm 
		% implementations. 
        %
        % Motivation: Migrate to a consistent coordinate definition and
        % fully vectorized computational engine. This implementation will
        % be fully backwards compatible with existing array definitions. It
        % will simply update the definition internally and expose new
        % methods to the user with the prefix "upd_"
        
        % Todo list
        % Consolidate to only FEKO element pattern (2D) format.
        function upd_redefineCoordinates(this)
            
           if this.upd_isvalid == 1
               disp('[AESA] Coordinates were previously redefined');
               return;
           end
                       
            % element position array shall be [N x 3]
            if(size(this.pos,2) ~= 3)
                this.pos = this.pos.';
                fprintf('Reshaped element position array to [%d x %d]\n',size(this.pos,1),size(this.pos,2));
            end
            
            % array normal shall be along (+) x dimension, y,z plane forms
            % standard planar array axis.
            if( any(this.pos(:,1) ~= 0 ) )
                this.pos = [this.pos(:,3) this.pos(:,1) this.pos(:,2)];
                fprintf('Updated frame definition to x(normal), y(up), z(left)\n');
            end
            
            vars = {'amp','ampaz','ampel','rng','phs','phs_steer','phs_sphcmp'};
            for k = 1 : length(vars)
                % Check "rng" vector
                var = vars{k};
                if( (size(eval(['this.' var]),1) ~= this.nElements) && ~isempty(eval(['this.' var])))
                    eval(['this.' var ' = this.' var '.'';']);
                    fprintf('Reshaped %s to [%d x %d]\n',var,size(eval(['this.' var]),1),size(eval(['this.' var]),2));
                end
                
            end
            
            this.upd_subArrayNumber = (1:this.nElements)';
            %this.upd_subArrayWgt = 1;
            this.upd_setSubArrayWeight(1); % Sets all sub-arrays to 1.
            this.upd_computeSubArrayCenters();
            
            this.upd_isvalid = 1; % set flag to allow updated calls
            
        end
        
        function voltageGain = upd_getGain(this,azRad,elRad,wgtMode,varargin)
            if( (nargin > 4) && strcmpi(varargin{1},'preserve_Az_El'))
                flag_optimize = 1;
            else
                flag_optimize = 0;
            end
            
            if(flag_optimize == 0) % full gain calculation (e.g. new az/el data)
                
                % Check if array is valid for updated functions
                if(this.upd_isvalid == 1)
                    
                    % If no weighting is requested, set to uniform
                    if(~exist('wgtMode'))
                        wgtMode = this.ENUM_MODE_UNIFORM;
                    end
                    
                    % Input expects a 1-D vector of azimuth & elevation points
                    nAz = length(azRad);
                    nEl = length(elRad);
                    this.upd_farField_staticAzEl = [nAz nEl];
                    
                    % Note: for indexing, M = # az pts * # el pts (nE*nA)
                    
                    % Az / El vectors must be [1 x (nE or nA)]
                    if(size(azRad,1) ~= 1)
                        azRad = azRad.';
                    end
                    if(size(elRad,1) ~= 1)
                        elRad = elRad.';
                    end
                    
                    % If single value passed in, meshgrid returns single
                    % value (e.g. it works regardless)
                    [azG,elG]=meshgrid(azRad,elRad);
                    
                    % convert requested az/el to unit vector
                    [ux,uy,uz]=sph2cart(azG,elG,1);
                    
                    % reform unit vector into [3 x (nAz * nEl)], or [3 x M]
                    R = [ ux(:)'
                          uy(:)'
                          uz(:)' ];
                    
                    % Compute phase map [N x 3] [3 x M]
                    P = this.twopi_ovr_lambda * this.pos * R ;
                    
                    if(this.upd_elementPatternPresent)
                        % Compute element pattern voltage gain as M x 1                
                        vGe = interp2( ...
                            this.upd_elementPatternAzG, ...
                            this.upd_elementPatternElG, ...
                            this.upd_elementPattern, ...
                            azG(:), ...
                            elG(:));
                        %apply peak gain from initialiation, or as manually
                        %defined
                        vGe = vGe.*10^(this.upd_elementGaindBV/10);
                    elseif(~isempty(this.upd_elementGaindBV))
                        vGe = 10^(this.upd_elementGaindBV/10);
                        vGe = (R(1,:) * vGe).';
                        vGe = vGe.'; %calculation for voltage gain expects column vector
                    else
                        error('No Element Pattern Defined. Either generate one or set fixed parameter upd_elementGaindBV');
                    end
                    
                    this.upd_farField_staticElement = vGe;
                    
                    % === Apply weighing to complex exponential and sum ===
                    wa = this.upd_getWeight(wgtMode);
                    
                    % Apply sub-array quadrature weights
                    % Note: this is a vector of ones unless specifically
                    % set by the user (e.g. has no effect by design)
                    wsub = this.upd_getSubArrayWeight();
                    
                    
                    % Steering Vector (set via upd_setSteer)
                    ws = exp(1i*this.phs_steer);
                    
                    % Build per element noise vector
                    this.upd_updateNoise();
                    
                    % Complete weighting vector [1 x N]
                    w = wa .* ws .* wsub.* this.upd_wn;
                    
                    % Compute gain (isotropic element) in voltage space [1 x M]
                     this.upd_farField_staticArray = exp(1i*P) ./ ...
                     sqrt(this.nElements); % KDS 12/12/2019 --- needs
                    % checked!!!
                    %this.upd_farField_staticArray = exp(1i*P) ./ this.nElements;
                    voltageGain = w.' * this.upd_farField_staticArray;
                    
                    % return the result in terms [El x Az] array format
                    % Note: if you just asked for a single az,el, you will
                    % get a single value back.
                    voltageGain = sqrt(this.efficiency) * ...
                        reshape(voltageGain .* vGe.',nEl,nAz) ;
                else
                    error('Updated array flag not set. Exiting.');
                end
            else % optimized (assumes already called and just updating weightings)
                
                % Steering Vector (set via upd_setSteer)
                ws = exp(1i*this.phs_steer);
                wa = this.upd_getWeight(wgtMode);
                
                % Build per element noise vector
                this.upd_updateNoise();
                
                wsub = this.upd_getSubArrayWeight();


                % Complete weighting vector [1 x N]
                w = wa .* ws .*  wsub.* this.upd_wn;
                
                % Static array already computed. Just leverage new
                % weightings.
                voltageGain = w.' * this.upd_farField_staticArray;
                    
                % return the result in terms [El x Az] array format
                % Note: if you just asked for a single az,el, you will
                % get a single value back.
                voltageGain = sqrt(this.efficiency) * ...
                    reshape(voltageGain .* this.upd_farField_staticElement.', ...
                    this.upd_farField_staticAzEl(2), ...
                    this.upd_farField_staticAzEl(1));
                
            end
        end
        
        function upd_setSteer(this,azRad,elRad)
            if(this.upd_isvalid == 1)
                this.phs_steer = this.upd_setSteerArray(azRad,elRad);
                %this.azSteer = azRad;
                %this.elSteer = elRad;
                % Build steering vector components
                %[uxs,uys,uzs]=sph2cart(azRad,elRad,1);

                % Form steering vector
                %rs = [uxs uys uzs]';

                % Compute phasing (rad) over array elements [N x 1]
                %this.phs_steer = -this.twopi_ovr_lambda * this.pos * rs;
                
                if(this.pShiftQuantization == 1)
                    this.needPhaseQuantized = 1;
                    this.upd_applyPhaseQuantization();
                end
                
                this.upd_ok2updateNoise = 1;
            else
                error('Updated array flag not set. Exiting.');
            end
        end
        
        function phs_steer_array = upd_setSteerArray(this,azRad,elRad,steerMode)
            % Purpose: Compute steering vectors (in rad) for array of az,el inputs
            %          steerMode = 'full' <default> or 'subarray'
            % NOTE: This function does not set any array internals
            % Each column represents a specific steering vector (rad) for
            % [A1E1,A1E2,...,A1EnEl,A2E1,A2,E2,...,AnAzEnEl]
            
            if(~exist('steerMode','var'))
                steerMode = 'full';
            end
            
            if(strcmpi(steerMode,'subarray'))
                % compute the subarray phase centers
                this.upd_computeSubArrayCenters();
            end

            % Input expects a 1-D vector of azimuth & elevation points
            nAz = length(azRad);
            nEl = length(elRad);

            % Note: for indexing, M = # az pts * # el pts (nE*nA)

            % Az / El vectors must be [1 x (nE or nA)]
            if(size(azRad,1) ~= 1)
                azRad = azRad.';
            end
            if(size(elRad,1) ~= 1)
                elRad = elRad.';
            end

             % If single value passed in, meshgrid returns single
            % value (e.g. it works regardless)
            [azG,elG]=meshgrid(azRad,elRad);

            % Build steering vector components
            [uxs,uys,uzs]=sph2cart(azG(:),elG(:),1);
            
            rs = [uxs uys uzs]';
            
            if(strcmpi(steerMode,'subarray'))
                phs_steer_array = -this.twopi_ovr_lambda * this.upd_subArrayPos * rs;
            else
                phs_steer_array = -this.twopi_ovr_lambda * this.pos * rs;
            end

        end
        
        function wa = upd_getWeight(this,wgtMode)
            % Determine amplitude weighting to use (file-based)
            switch(wgtMode)
                case this.ENUM_MODE_UNIFORM
                    wa = ones(this.nElements,1);
                case this.ENUM_MODE_AMPFILE
                    wa = this.amp;
                case this.ENUM_MODE_DAZ_AMPFILE
                    wa = this.ampaz;
                case this.ENUM_MODE_DEL_AMPFILE
                    wa = this.ampel;
                case this.ENUM_MODE_ABF
                    wa = this.abf_term;
                otherwise
                    error('Unknown weigting mode passed in. Exiting.');
            end
            
            if( isempty(wa) )
                error(sprintf('Weighting %d is empty! Exiting!',wgtMode));
            end 
        end
        
        function upd_setSubArrayWeight(this,wgt)
            nSubArrays = length(unique(this.upd_subArrayNumber));
            if(~exist('wgt','var'))
                wgt = 1;
                fprintf('Weights reset to 1.0\n');
            end
            if( (wgt(1) == 1) && (length(wgt) == 1) )
                this.upd_subArrayWgt = ones(nSubArrays,1);
                
                
            else
                if(size(wgt,1) < size(wgt,2))
                    wgt = wgt.';
                end
                
                if(length(wgt) ~= nSubArrays)
                    error('There are %d sub-arrays defined, and only %d weights passed in!',nSubArrays,length(wgt));
                end
                this.upd_subArrayWgt = wgt;
            end
            this.upd_genSubArrayMatrix();
        end
        
        function upd_computeSubArrayCenters(this)
            this.upd_subArrayPos = [];
            iSubArrays = unique(this.upd_subArrayNumber);
            nSubArrays = length(iSubArrays);
            for k = 1 : nSubArrays
                p=mean(this.pos(this.upd_subArrayNumber==iSubArrays(k),:),1);
                this.upd_subArrayPos(k,:) = p;
            end
        end
        
        function wa = upd_getSubArrayWeight(this,varargin)
            % Returns current sub-array weights
            % Inputs are in name,value pairs
            % weight_vector = <obj>.upd_getSubArrayWeight(<mode>,<subarray #s>)
            %  <mode>: 'mode', 'fullarray' or 'subarray'
            %        -full array option will return sub array weights as
            %        they apply to full array 
            %        -sub array option will return just the sub-array
            %        weights
            %        ** default is 'fullarray'
            %
            %  <subarray #s>: 'subarrays', [sub array numbers of interest]
            %        ** default is all subarrays
            
            nvarargin = length(varargin);
            
            % set up defaults
            amode = 1; % output full repeated array weighting of s.a. weights
            saNum = unique(this.upd_subArrayNumber); % all sub arrays
            
            for kv = 1 : 2 : nvarargin
                switch lower(varargin{kv})
                    case lower('mode')
                        switch lower(varargin{kv+1})
                            case lower('subarray')
                                amode = 0;
                            case lower('fullarray')
                                amode = 1;
                            otherwise
                                error('Unknown array mode passed in (%s)\n',varargin{kv+1});
                        end
                    case lower('subarrays')
                        saNum = varargin{kv+1};
                    otherwise
                        error('Unknown option: %s\n',varargin{kv});
                end
            end
            
            % make sure the maximum requested sub-array weight is defined
            maxSa = max(saNum);
            if(length(this.upd_subArrayWgt) < maxSa)
                error('Subarray %d is requested but largest defined weight is %d!\n',maxSa,length(this.upd_subArrayWgt));
            end
            
            
            if(amode == 1)
            
            
            % Ensure sub-array mapping matrix is updated
            this.upd_genSubArrayMatrix();
            
            S = this.upd_subArrayMatrix; % matrix is [Ne x Msub]
            
            % 
            Sreq =0*S;
            
            for k = 1 : length(saNum)
                try
                    Sreq(:,saNum(k)) = S(:,saNum(k));
                catch
                    if(~any(this.upd_subArrayNumber==saNum(k)))
                        error(sprintf('Sub-Array %d not found! Exiting!',saNum(k)));
                    else
                        error(lasterr);
                    end
                end
            end
            wa = Sreq * this.upd_subArrayWgt;
            else % user is requesting just sub-array weights
                wa = this.upd_subArrayWgt(saNum);
            end
            
        end
        
        function upd_updateNoise(this)
            % ============== Noise Processing =====================
            % Updates noise (if enabled) in the object (upd_wn)
            % Only update noise if this flag is set high. User code must
            % reset this flag high to generate new noise.
            if(this.upd_ok2updateNoise)
                
                % If noise is enabled, it applies to the N element
                % space as a weighting
                wa_noise = 1;
                wp_noise = 0;
                
                if(this.aVar ~= 0)
                    wa_noise = 1 + this.aVar * randn(this.nElements,1);
                end
                
                if(this.pVar ~= 0)
                    wp_noise = 0 + this.pVar * randn(this.nElements,1);
                end
                
                this.upd_wn = wa_noise .* exp(1i*wp_noise);
                this.upd_ok2updateNoise = 0;
                % ============== Noise Processing =====================
            end
        end
        
        function upd_scaleWeights(this)
            if this.upd_weightsScaled
                disp('[AESA] Weights were previously scaled');
                return;
            end
            
            p_uniform = this.computeCumulativeNoise(this.ENUM_MODE_UNIFORM);
            if(~isempty(this.amp))
                p_amp = this.computeCumulativeNoise(this.ENUM_MODE_AMPFILE);
                this.amp = this.amp .* sqrt(p_uniform/p_amp);
            end
            if(~isempty(this.ampaz))
                p_ampaz = this.computeCumulativeNoise(this.ENUM_MODE_DAZ_AMPFILE);
                this.ampaz = this.ampaz .* sqrt(p_uniform/p_ampaz);
            end
            if(~isempty(this.ampel))
                p_ampel = this.computeCumulativeNoise(this.ENUM_MODE_DEL_AMPFILE);
                this.ampel= this.ampel .* sqrt(p_uniform/p_ampel);
            end
            
            this.upd_weightsScaled = 1;
        end
        
        function [hp,G] = upd_plot(this,azPts,elPts,wgtModes,varargin) % | wgt_str
            if(nargin == 1)
                azPts = linspace(-pi/2,pi/2,400);
                elPts = linspace(-pi/2,pi/2,400);
                wgtModes = this.ENUM_MODE_UNIFORM;
                varargin{1} = 'Default Uniform';
                warning('plotting using defaults!');
                fprintf('[hp,G] = upd_plot(this,azPts(r),elPts(r),wgtModes,varargin)\n');
            end
            
            azCutFlag = 0;
            elCutFlag = 0;
            fullFlag = 0;
            
            % wgtModes can be a vector of weight modes, and will lead to
            % multiple plots.
            
            % determine if az cut, el cut, or full pattern
            if(length(azPts) == 1)
                elCutFlag = 1;
                title_str = sprintf('%s Elevation Cut @ %4.1f\\circ Azimuth (dBi)',this.name,180/pi*azPts);
            elseif(length(elPts) == 1)
                azCutFlag = 1;
                title_str = sprintf('%s Azimuth Cut @ %4.1f\\circ Elevation (dBi)',this.name,180/pi*elPts);
            else
                fullFlag = 1;
                title_str = sprintf('%s Far Field Pattern (dBi)',this.name);
            end
            
            for kwgt = 1 : length(wgtModes)
                wgtMode = wgtModes(kwgt);
            
            if(nargin == 4) % no weighting string passed in
                % Determine weighting string
                switch(wgtMode)
                    case this.ENUM_MODE_UNIFORM
                        wgt_str = 'Uniform Weighting';
                    case this.ENUM_MODE_AMPFILE
                        wgt_str = 'Amplitude File';
                    case this.ENUM_MODE_DAZ_AMPFILE
                        wgt_str = 'Delta Azimuth File';
                    case this.ENUM_MODE_DEL_AMPFILE
                        wgt_str = 'Delta Elevation File';
                    case this.ENUM_MODE_ABF
                        wgt_str = 'Adaptive';
                    otherwise
                        error('Unknown weighting!');
                end
                if(length(wgtModes) > 1)
                    warning('If you are using multi-plot, please pass in {a,b,...} legend string');
                end
            else
                wgt_str = varargin{1};
            end
            this.upd_wgt_name = wgt_str;
            
            
            
            % request gain data
            vGain = this.upd_getGain(azPts,elPts,wgtMode);
            G = 20*log10(abs(vGain));
            Gbounds = [max(G(:))-70  max(G(:))+5];
            azVec = azPts * 180/pi;
            elVec = elPts * 180/pi;
            if(fullFlag)
                
                if( length(wgtModes) > 1 )
                    figure;
                end
                
                hi=imagesc(elVec, azVec, G);
                set(gca,'YDir','normal');
                xlabel('Azimuth (deg)');
                ylabel('Elevation (deg)');
                title({title_str,wgt_str});
                colorbar;
                caxis(Gbounds);
                xlim([min(azVec) max(azVec)]);
                ylim([min(elVec) max(elVec)]);
                f1pos = get(gcf,'Position');
                if( length(wgtModes) == 1)
                [px,py,pz] = sphere(50);
                
                figure;
                sPattern = surface(py, px ,(pz));  
                sPattern.FaceColor = 'texturemap';        % set color to texture mapping
                sPattern.EdgeColor = 'none';              % remove surface edge color
                sPattern.CData = [get(hi,'CData') 0*get(hi,'CData')-999];
                

                set(gca,'DataAspectRatio',[1 1 1],'Color','k', ...
                    'XColor','w','YColor','w','ZColor','w','YDir','reverse');
                ylabel('U Space');
                zlabel('V Space');
                hc=colorbar;
                set(hc,'Color','w')
                set(gcf,'Color','k');
                title({title_str,wgt_str},'Color','w');
                caxis(Gbounds);
                view(-90,0);
                this.upd_handle_2d_plots = [hi sPattern];
                f2pos = get(gcf,'Position');
                f2pos(1) = f1pos(1)+f1pos(3);
                set(gcf,'Position',f2pos);
                hp = hi;
                
                if(isempty(get(hi,'ButtonDownFcn')))
                    set(hi,'ButtonDownFcn',{@this.upd_update2D,wgtMode,azPts,elPts},'UserData',wgtMode);
                    c = uicontextmenu;
                    set(c,'Parent',get(get(hi,'Parent'),'Parent'));

                    % Assign the uicontextmenu to the plot line
                    hi.UIContextMenu = c;

                    % Create child menu items for the uicontextmenu
                    m1 = uimenu(c,'Label','Uniform','Callback', {@this.upd_update2D,0,azPts,elPts});
                    m2 = uimenu(c,'Label','Weighted','Callback',{@this.upd_update2D,1,azPts,elPts});
                    m3 = uimenu(c,'Label','Delta Az','Callback',{@this.upd_update2D,2,azPts,elPts});
                    m4 = uimenu(c,'Label','Delta El','Callback',{@this.upd_update2D,3,azPts,elPts});
                    m5 = uimenu(c,'Label','Adaptive','Callback',{@this.upd_update2D,4,azPts,elPts});

                end
                
                end
                
                
            elseif(azCutFlag)
                if(nargout < 2)
                    hp=plot(azVec,G);hold on;
                    grid on;
                    xlabel('Azimuth (deg)');
                    if(length(wgtModes) == 1)
                        title({title_str,wgt_str});
                    elseif(kwgt == length(wgtModes))
                        title(title_str);
                        legend(wgt_str);
                    end
                    ylim(Gbounds);
                else
                    datOut = G;
                    hp = [];
                end
            elseif(elCutFlag)
                if(nargout < 2)
                    hp=plot(elVec,G);hold on;
                    grid on;
                    xlabel('Elevation (deg)');
                    if(length(wgtModes) == 1)
                        title({title_str,wgt_str});
                    elseif(kwgt == length(wgtModes))
                        title(title_str);
                        legend(wgt_str);
                    end
                    ylim(Gbounds);
                else
                    datOut = G;
                    hp = [];
                end
            end
            
            end
            
        end
        
        function upd_update2D(this,varargin)
            
            if(nargin > 1) % got here via callback
                
                % if event is hit and button is 1, re-steer + re-draw
                % if event is action re-draw only
                redraw = 0;
                if( strcmp(varargin{2}.EventName,'Hit') )
                    % Only re-steer if primary clicked
                    if(varargin{2}.Button == 1)
                        if( ...
                            sum(((real(this.upd_subArrayWgt) == 1) & ...
                            (imag(this.upd_subArrayWgt) == 0))) == 0)
                        errordlg('Adaptive weights have been computed. No steering allowed via click.','No can do my friend...');
                        return;
                        end
                        cp = get(gca,'CurrentPoint');
                        azSteer_click = cp(1,1)*pi/180;
                        elSteer_click = cp(1,2)*pi/180;
                        this.upd_setSteer(azSteer_click,elSteer_click);
                        redraw=1;
                    end
                elseif( strcmp(varargin{2}.EventName,'Action') )
                    if(isvalid(this.upd_handle_2d_plots(1)))
                        set(this.upd_handle_2d_plots(1),'UserData',varargin{3});
                        switch(varargin{3})
                            case 0
                                this.upd_wgt_name = 'Uniform';
                            case 1
                                this.upd_wgt_name = 'Amp File';
                            case 2
                                this.upd_wgt_name = 'DAz File';
                            case 3
                                this.upd_wgt_name = 'DEl File';
                            otherwise
                                this.upd_wgt_name = 'Adaptive';
                        end
                    else
                        warning('2D Plot Handle Invalid. Skipping');
                    end
                    redraw=1;
                end
            else % called externally
                % assume external script has updated everything.
                redraw=1;
            end
            
            if(redraw)
                
                if(isvalid(this.upd_handle_2d_plots(1)))
                    wgtMode = get(this.upd_handle_2d_plots(1),'UserData');
                end
                
                vGain = this.upd_getGain([],[],wgtMode,'preserve_Az_El');
                G = 20*log10(abs(vGain));
                
                if(isvalid(this.upd_handle_2d_plots(1)))
                    set(this.upd_handle_2d_plots(1),'CData',G);
                else
                    warning('2D Plot Handle Invalid. Skipping');
                end
                
                if(isvalid(this.upd_handle_2d_plots(2)))
                    set(this.upd_handle_2d_plots(2),'CData',[G 0*G-999]);
                else
                    warning('Surface Plot Handle Invalid. Skipping.')
                end
                
            end
        end
        
        function upd_genSubArrayMatrix(this)
            % find out how many subarrays there are
            iSubs = unique(this.upd_subArrayNumber);
            nSubs = length(iSubs);
            
            % dimension sub array matrix to [Ne x Nsubs]
            this.upd_subArrayMatrix = zeros(this.nElements,nSubs);
            
            for k = 1 : nSubs
                % set ones in the sub array matrix appropriatatly 
                this.upd_subArrayMatrix( ...
                    this.upd_subArrayNumber == iSubs(k), ...
                    k ) = 1;
            end
            
        end
        
        function ue = upd_genIQ_Channels(this,u,varargin)
            % u         input time domain waveform [1 x # samples]
            % fromAz    azimuth (rad) of signal w.r.t. array normal
            % fromEl    elevation (rad) of signal w.r.t. array normal
            % subArray  [optional] output only elements who are part of sub
            %           array # <subArray>. All elements default to group 1
            % applySteering     [varargin] 1 = apply phs_steer to signals,
            %                   0 = direct elemental sampling
            % applyWeight       maps to enumerated weighting
            % sumOutput         [varargin] 1 = apply 
            
            nvarargin = length(varargin);
            
            % set up defaults
            fromAz = 0;
            fromEl = 0;
            subArray = unique(this.upd_subArrayNumber);
            applySteering = 0;
            applyWeighting = 0;
            sumOutput = 0;
            
            for kv = 1 : 2 : nvarargin
                switch lower(varargin{kv})
                    case lower('fromAz')
                        fromAz = varargin{kv+1};
                    case lower('fromEl')
                        fromEl = varargin{kv+1};
                    case lower('subArray')
                        subArray = varargin{kv+1};
                    case lower('applySteering')
                        applySteering = varargin{kv+1};
                    case lower('applyWeighting')
                        applyWeighting = varargin{kv+1};
                    case lower('sumOutput')
                        sumOutput = varargin{kv+1};
                    otherwise
                        error(sprintf('Unknown option: %s\n',varargin{kv}));
                end
            end
            
            
            % Ensure time domain data is [1 x temporal samples]
            if(size(u,1) ~= 1)
                u = u.';
            end
            
            % Angle-of-arrival unit vector components
            [uxs,uys,uzs]=sph2cart(fromAz,fromEl,1);

            % Form angle-of-arrival unit vector
            rs = [uxs uys uzs]';

            % Compute phasing (rad) over array elements [N x 1]
            signalPhasing = this.twopi_ovr_lambda * this.pos * rs;
            
            % Apply element pattern
            if(this.upd_elementPatternPresent)
                % Compute element pattern voltage gain as M x 1                
                vGe = interp2( ...
                    this.upd_elementPatternAzG, ...
                    this.upd_elementPatternElG, ...
                    this.upd_elementPattern, ...
                    fromAz, ...
                    fromEl);
            elseif(~isempty(this.upd_elementGaindBV))
                vGe = 10^(this.upd_elementGaindBV/10);
            else
                error('No Element Pattern Defined. Either generate one or set fixed parameter upd_elementGaindBV');
            end
            
            if(applySteering)
                % Apply element gain and phasing to array of signals [N x TD]
                ws = exp(1i*this.phs_steer);
            else
                ws = ones(this.nElements,1);
            end
            
            % Apply element gain and phasing to array of signals [N x TD]
            wa = this.upd_getWeight(applyWeighting);
            
            % Get sub-array weight
            wsa = this.upd_getSubArrayWeight('mode','fullarray','subarrays',subArray);
            
            %   [1x1] *    [N x 1] * [1 x TD]
            ue_full = vGe .* (wa .* ws .* wsa .* exp(1i*signalPhasing)) * u;
            
            % there will be as many channels of output as there are
            % sub-arrays [Nsub x TD]
            ue = this.upd_subArrayMatrix(:,subArray)' * ue_full;
            
            % If user requests after the summing junction
            if(sumOutput)
                ue = sum(ue,1);
            end
            
            
            
            
                    
        end
        
        function Rest = upd_estCorrMtx(this,uArray,nSamples,subArray)
            % uArray - time domain signals per element [N x # time domain
            %          samples
            % nSamples - # of samples to use for correlation matrix
            %          estimation (empty uses all of them)
            % subArray - [optional] output only elements who are part of sub
            %            array # <subArray>. All elements default to group
            %            are in their own sub-array by default.
            
            % Apply subarray filtering if asked
            if(exist('subArray'))
                idx = find(this.upd_subArrayNumber == subArray);
                if(isempty(idx))
                    error(sprintf('Sub-Array %d not found. (<obj>.upd_subArrayNumber)',subArray));
                end
            else
                idx = 1 : size(uArray,1);
            end
            
            
            if(isempty(nSamples))
                nSamples = length(uArray(1,:));
            end
            
            Rest = uArray(idx,1:nSamples)*uArray(idx,1:nSamples)' ...
                   ./ nSamples;
            
        end
        
        function Dinv = diagInv(this,D,N)
            % This function takes the diagonal elements of a matrix, inverts them element-wise
            % and then returns them as a new diagonal matrix
            Dinv = diag(1.0 ./ diag(D));
            
            if(exist('N','var'))
                for k = (N+1):length(D)
                    Dinv(k,k) = 0;
                end
            end

        end


% Direct Az,El Method
%%
% [ag,eg]=meshgrid((-90:.5:90)*pi/180,(-90:.5:90)*pi/180);
% ag = single(ag);
% eg = single(eg);
% [uz,ux,uy] = sph2cart(ag,eg,1);
% uvec = [ux(:) uy(:) uz(:)];
%
% % compute terms
% ovec = ones(length(uvec),1);
% this_a = 1 + p2.aVar * randn(1,p2.nElements,p2.precisionClass);
% this_p = 0+ p2.pVar * randn(1,p2.nElements,p2.precisionClass);
% p2.phs = (uvec * p2.pos) * p2.twopi_ovr_lambda + ovec * (this_p - p2.phs_steer);
%
% val = exp(1i*p2.phs);
% csum = (ovec * this_a              ) .* val;
% csum = sum(csum,2);
% this_elFac = 1.0; % todo fix this for large az/el calc.
% this_elFac=sqrt(abs(uz(:)));
% gain = csum .* p2.neScale .* p2.efficiency .* p2.vPeak .* this_elFac;
% g3 = reshape(gain.',size(ag,1),size(ag,2));
%
% figure;
% imagesc(20*log10(abs(g2)));
% colorbar;
%
% figure;
% imagesc(20*log10(abs(g3)));
% colorbar;

function importFEKOElement(this,fileName,options)
%Author: Joel Schopis Senior Consultant Booz Allen Hamilton
%This function reads a FEKO patterns and assigns it into the
%element pattern.  Pattern is converted to Azimuth over Elevtion
%polarization.  Elevation is considered vertical and Azimuth
%horizontal.
%Variables
%   this: The current array object
%   fileName: FEKO file from which to read the radiation pattern,
%       radiation pattern must be a sphere in range.
%   options: struct that contains options for pattern:
%       options.copol: string that deterimines which polarization will
%       be copol, is case insensitive default is horizontal
%       possible values:
%           Vertical: ('el', 'vertical', 'v')
%           Horizontal: ('az','horizontal', 'h')
%           RHCP: ('rhcp','right','r')
%           LHCP: ('lhcp','left', 'l')
%       options.crosspol: string that deterimines which polarization will
%       be copol, is case insensitive defualt is vertical
%       possible values are same as copol variable
%       options.increment: degrees, increment over element
%       radiation pattern, default is 1 degree.

copol = 'horizontal'; %Copolarization
crosspol = 'vertical'; %Crosspolarization
patIncr = 1; %Element pattern increment
rotM = []; %Rotation matrix for pattern, currently not in use

%Set options
if(nargin > 2 && ~isempty(options))
    %Set copolarization option
    if(isfield(options,'copol'))
        switch(lower(options.copol))
            case {'el','vertical','v'}
                copol = 'vertical';
            case {'az','horizontal','h'}
                copol = 'horizontal';
            case {'rhcp','right','r'}
                copol = 'rhcp';
            case {'lhcp','left','l'}
                copol = 'lhcp';
            otherwise
                warning(sprintf('Given copolarization ''%s'' is not valid type, using horizontal polarization.',options.copol));
                copol = 'horizotal';
        end
    end
    
    %Set crosspolarization option
    if(isfield(options,'crosspol'))
        switch(lower(options.copol))
            case 'el'
            case 'vertical'
            case 'v'
                crosspol = 'vertical';
            case 'az'
            case 'horizontal'
            case 'H'
                crosspol = 'horizontal';
            otherwise
                warning(sprintf('Given crossopolarization ''%s'' is not valid type, using vertical polarization.',options.crosspol));
                crosspol = 'vertical';
        end
    end
    
    %Set element increment
    if(isfield(options,'increment'))
        patIncr = options.increment;
    end
    
    if(isfield(options,'rotate') && numel(options.rotate))
        rotM = eye(3,3);
        rotMrev = eye(3,3);
        
        for idx = 1:size(options.rotate,1)
            rotConfig = options.rotate(idx,:);
            
            switch(rotConfig(1))
                case 1
                    %Rotate around X axis
                    rot = [ 1 0 0;
                        0 cos(rotConfig(2)) -sin(rotConfig(2));
                        0 sin(rotConfig(2))  cos(rotConfig(2));];
                    
                    rotRev = [ 1 0 0;
                        0 cos(-rotConfig(2)) -sin(-rotConfig(2));
                        0 sin(-rotConfig(2))  cos(-rotConfig(2));];
                case 2
                    %Rotate around Y axis
                    rot = [ cos(rotConfig(2)) 0 sin(rotConfig(2));
                        0 1            0;
                        -sin(rotConfig(2)) 0 cos(rotConfig(2));];
                    
                    rotRev = [ cos(-rotConfig(2)) 0 sin(-rotConfig(2));
                        0 1            0;
                        -sin(-rotConfig(2)) 0 cos(-rotConfig(2));];
                case 3
                    %Rotate around Z axis
                    rot = [ cos(rotConfig(2)) -sin(rotConfig(2)) 0;
                        sin(rotConfig(2))  cos(rotConfig(2)) 0;
                        0  0  1;];
                    
                    rotRev = [ cos(-rotConfig(2)) -sin(-rotConfig(2)) 0;
                        sin(-rotConfig(2))  cos(-rotConfig(2)) 0;
                        0  0  1;];
                otherwise
                    error('Rotation type not found');
            end
            
            rotM = rotM*rot;
            rotMrev = rotRev*rotMrev;
        end
    end
end

%Open FEKO file
fid = fopen(fileName,'r');
if(~fid)
    warning('Unable to open file. Cowardly refusing to set element pattern.');
    return;
end

feko_data = zeros(100000,9); %Intialize feko data to increase run time
count = 0; %Count the number of FEKO rows with data that are read.
exportRead = false;
exportType = 'Unkown';

%Loop through file to get feko data
while(~feof(fid))
    %Get export data
    data = fgetl(fid);
    
    %Check export type, too see if theta or phi is incremented
    %first
    if(~exportRead)
        %Get what program generated far field
        feko_row = sscanf(data,'** File exported by %s');
        
        if(~isempty(feko_row))
            switch(upper(feko_row))
                case 'FEKO'
                    %Farfield was generated by FEKO kernal
                    exportRead = true;
                    exportType = 'FEKO';
                case 'POSTFEKO'
                    %Farfield was generated by POSTFEKO
                    exportRead = true;
                    exportType = 'POSTFEKO';
                otherwise
                    %Farfield was generated by POSTFEKO
                    exportRead = true;
                    exportType = 'FEKO';
                    warning('Could not determine how far field file was generated.  Assuming FEKO kernal.');
            end
        end
    else
        %Search for FEKO row that follows the appropriate format
        %skip all other rows.
        feko_row = sscanf(data,'%f');
        if ~isempty(feko_row)
            count = count + 1;
            
            %FEKO data variable is full double feko variable size,
            %this will increase speed of loop.
            if(count > length(feko_data))
                temp_data = zeros(length(feko_data)*2,9);
                temp_data(1:(count - 1),:) = feko_data(1:(count - 1),:);
                feko_data = temp_data;
            end
            
            feko_data(count,:) = feko_row(:);
        end
    end
end
fclose(fid);

%Remove excessive rows in FEKO data variable
feko_data = feko_data(1:count,:);

%Set theta range
THETA = unique(feko_data(:,1))*pi/180;
%set phi range
PHI = unique(feko_data(:,2))*pi/180;

%Put FEKO data into complex format and in PHI/THETA
%polarization
Theta_data = feko_data(:,3) + 1j*feko_data(:,4);
Phi_data = feko_data(:,5) + 1j*feko_data(:,6);

Theta_data(Theta_data == 0) = 1e-10;
Phi_data(Phi_data == 0) = 1e-10;

switch(exportType)
    case 'FEKO'
        Theta_data = reshape(Theta_data,[length(THETA) length(PHI)]);
        Phi_data = reshape(Phi_data,[length(THETA) length(PHI)]);
    case 'POSTFEKO'
        Theta_data = reshape(Theta_data,[length(PHI) length(THETA)])';
        Phi_data = reshape(Phi_data,[length(PHI) length(THETA)])';
    otherwise
        Theta_data = reshape(Theta_data,[length(THETA) length(PHI)]);
        Phi_data = reshape(Phi_data,[length(THETA) length(PHI)]);
end

%Grid format for Phi/Theta coordinates
[PHI,THETA] = meshgrid(PHI,THETA);

%             keyboard
%Azimuth and Elevation coordinates for element pattern
Az = (-180:patIncr:180)*pi/180;
El = (-90:patIncr:90)*pi/180;

%Turn element azimuth and elevation vectors into grids
[azGrid elGrid] = meshgrid(Az,El);

%Convert the element Az/El coordinates to PHI/THETA
THETAae = acos(cos(elGrid).*cos(azGrid));
PHIae = atan2(sin(elGrid),cos(elGrid).*sin(azGrid));

%Obtain the element PHI/THETA radiation pattern from FEKO
%patterns
THETAaePol = interp2(PHI,THETA,Theta_data,PHIae,THETAae);
PHIaePol = interp2(PHI,THETA,Phi_data,PHIae,THETAae);

%Convert element PHI/THETA radiation pattern into cartesian
%polarization
UaePol = -PHIaePol.*sin(PHIae) + THETAaePol.*cos(PHIae).*cos(THETAae);
VaePol = PHIaePol.*cos(PHIae) + THETAaePol.*sin(PHIae).*cos(THETAae);
WaePol = -THETAaePol.*sin(THETAae);

%Rotate FEKO pattern, not yet functional
if(~isempty(rotM))
    %Cartesian coordinates for FEKO patterns
    U = cos(elGrid).*sin(azGrid);
    V = sin(elGrid);
    W = cos(elGrid).*cos(azGrid);
    cartCoord = [U(:) V(:) W(:)];
    
    clear U V W;
    
    %Rotate coordinates
    rotCartCoord = rotM*cartCoord';
    clear cartCoord;
    rotU = reshape(rotCartCoord(1,:),size(elGrid));
    rotV = reshape(rotCartCoord(2,:),size(elGrid));
    rotW = reshape(rotCartCoord(3,:),size(elGrid));
    clear rotCartCoord;
    
    %Convert back to sphere and put back in grid form
    rotAz = atan2(rotU,rotW);
    rotEl = asin(rotV);
    
    clear rotU rotV rotW;
    
    %Get electric field vectors for rotation
    UPol = interp2(Az,El,UaePol,rotAz,rotEl);
    VPol = interp2(Az,El,VaePol,rotAz,rotEl);
    WPol = interp2(Az,El,WaePol,rotAz,rotEl);
    
    %Rotate electric field vectors
    cartPol = [UPol(:) VPol(:) WPol(:)];
    clear UPol VPol WPol;
    
    rotCartPol = rotMrev*cartPol';
    UaePol = reshape(rotCartPol(1,:),size(elGrid));
    VaePol = reshape(rotCartPol(2,:),size(elGrid));
    WaePol = reshape(rotCartPol(3,:),size(elGrid));
    clear rotCartPol;
end

%Store cartesian polarization parts for possible later use
this.elementPolCart.U = UaePol;
this.elementPolCart.V = VaePol;
this.elementPolCart.W = WaePol;

%Convert element cartesian radiation pattern into AZ/EL
%polarization.  This is the final polarization for the element
%pattern.
AzPol = UaePol.*cos(azGrid) - WaePol.*sin(azGrid);
ElPol = -UaePol.*sin(elGrid).*sin(azGrid) + ...
    VaePol.*cos(elGrid) - ...
    WaePol.*sin(elGrid).*cos(azGrid);

%Set element points
this.elemAz = azGrid;
this.elemEl = elGrid;

%Calculate Circular polarization
RhcPol = (AzPol + 1i*ElPol)/sqrt(2);
LhcPol = (AzPol - 1i*ElPol)/sqrt(2);

%Normalize pattern
AzPol = AzPol./max(AzPol(:));
ElPol = ElPol./max(ElPol(:));
RhcPol = RhcPol./max(RhcPol(:));
LhcPol = LhcPol./max(LhcPol(:));

%Assign copolar
if strcmpi('rhcp',copol)
    this.elementPattern = RhcPol;
elseif strcmpi('lhcp',copol)
    this.elementPattern = LhcPol;
elseif strcmpi('vertical',copol)
    this.elementPattern = ElPol;
else
    this.elementPattern = AzPol;
end

%Assign crosspol
if(isfield(this.Settings,'crossPol') && this.Settings.crossPol)
    if strcmpi('rhcp',crosspol)
        this.elementPatternCross = RhcPol;
    elseif strcmpi('lhcp',crosspol)
        this.elementPatternCross = LhcPol;
    elseif strcmp('vertical',crosspol)
        this.elementPatternCross = ElPol;
    else
        this.elementPatternCross = AzPol;
    end
end
end

    end
end


