classdef array < handle
    properties
        ENUM_MODE_UNIFORM = 0;
        ENUM_MODE_AMPFILE = 1;
        ENUM_MODE_DAZ_AMPFILE = 2;
        ENUM_MODE_DEL_AMPFILE = 3;
        ENUM_MODE_ABF = 4;
        pos
        elem_area
        elemfac_gain
        elemfac_ang
        elemfac_az
        elemfac_el
        amp,ampel,ampaz
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
        elementPattern
        vPeak,aVar,pVar,azSteer,elSteer,azSteerThresh,elSteerThresh
        adcSample
        pShiftQuantization,needPhaseQuantized
        pMin,pMax,pBits,pStep
        
        acBodyToAntenna = eye(3); % assume antenna straight out nose
        
        % matlab specific
        calcMethod = 2; % 2 = vectorized, 1 = exact FORTRAN copy
        gpuMethod = 0;
        precisionClass = 'single';
        name = '';
        plot_marker_size = 32;
        plot_2d_res_deg = 1.0;
        
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
        
        function generateIQData(this,steerAzRad,steerElRad,fromAzRad,fromElRad,refIQ)
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
            if(isempty(this.wavelength))
                error('wavelength field not defined!!');
            end
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
                ampfile_az,ampfile_el,elfacfile)
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
                rawdata = fscanf(fid,'%f');
                fclose(fid);
                this.amp = rawdata(2:end).';
            end
            if( ~isempty(ampfile_az) )
                fid = fopen(ampfile_az,'r');
                rawdata = fscanf(fid,'%f');
                fclose(fid);
                this.ampaz = rawdata(2:end).';
            end
            if( ~isempty(ampfile_el) )
                fid = fopen(ampfile_el,'r');
                rawdata = fscanf(fid,'%f');
                fclose(fid);
                this.ampel = rawdata(2:end).';
            end
            if( ~isempty(elfacfile) )
                % Reference syntax:
                % this_elFac = interp2(this.elemfac_az,this.elemfac_el,this.elemfac_gain,az,el);
                
                fid = fopen(elfacfile,'r');
                % Assumed format
                % nAz (integer)
                % nEl (integer)
                % el pts (complex, no imaginary)
                % az pt cgain,cgain, etc
                nAz = fscanf(fid,'%d',1);
                nEl = fscanf(fid,'%d',1);
                rawData = zeros(2,nAz*nEl + nAz + nEl);
                fgetl(fid);
                k = 1;
                while(~feof(fid))
                    s = fgetl(fid);
                    rawData(:,k) = sscanf(s,'%f,%f');
                    k = k + 1;
                end
                fclose(fid);

                elPts = rawData(1,1:nEl);
                azGnPts = rawData(1,(nEl+1):end) + 1i * rawData(2,(nEl+1):end);

                %%
                G = reshape(azGnPts,nEl+1,nAz).';
                azPts = G(:,1);
                gPts = G(:,2:end);
                
                [azG,elG]=meshgrid(azPts,elPts);
                this.elemfac_az = azG;
                this.elemfac_el = elG;
                this.elemfac_gain = gPts.';
                this.elemFacPresent = 1;
                
                this.vPeak = max(abs(gPts(:))) * sqrt(this.efficiency) * this.nElements;
                
            end
            
        end % init stap array
        
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
        
        function preAllocate(this)
            this.pos = zeros(3,this.nElements,this.precisionClass);
            this.amp = zeros(1,this.nElements,this.precisionClass);
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
            % if you pass in a 2nd argument, it will color code
            % the elements.
            %hp = plot3(this.pos(1,:),this.pos(2,:),this.pos(3,:),'.');
            if(nargin == 1)
                hs=scatter3( ...
                this.pos(1,:), ...
                this.pos(2,:), ...
                this.pos(3,:), ...
                this.plot_marker_size);
                set(hs,'Marker','sq','MarkerFaceColor','r');
            else
                hs=scatter3( ...
                this.pos(1,:), ...
                this.pos(2,:), ...
                this.pos(3,:), ...
                this.plot_marker_size, ...
                varargin{1},'filled');
                set(hs,'Marker','sq');
                colorbar;
            end
            view(0,90);
            xlabel('X Direction (m)');
            ylabel('Y Direction (m)');
            zlabel('Z Direction (m)');
            title('Array Geometry');
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
        
        function gain = getGain(this,azInertial,elInertial,mode)
             
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
                disp('phase quant not implemented');
            end
            
            % compute element factor for this far-field location
            if(this.elemFacPresent == 1)
                %this_elFac = GetInterp_D(this.elementPattern,incAngle);
%                 this_elFac = interp1(this.elemfac_ang,this.elemfac_gain,incAngle);
                this_elFac = interp2(this.elemfac_az,this.elemfac_el,this.elemfac_gain,az,el);
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
                this_a = this_a .* (~this.element_off);
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
                elseif(mode == this.ENUM_MODE_ABF)
                    if(isempty(this.abf_term))
                        error('you must set member property: abf_term');
                    end
                    
                    %%%%%%%%%%val = exp(1i*this.phs) .* this.abf_term; % convert to complex weights and
                    % apply adapative beamformer weights
                    % csum = (this_a .* this.abf_term) .* val .* exp(1i*this.phs_steer);
                    csum = onecolvec*(this_a .* this.abf_term) .* val;
                end
                %csum = sum(csum); % JAH - original
                csum = sum(csum,2); % JAH - vectorized inputs, 2nd arg sum over N-elements dimensions (sum each row)
                csum = reshape(csum,1,[]); % JAH - vectorized inputs, put back into ROW VECTOR
            end
            
            % jcl(-)gain = csum * this.neScale * this.efficiency * this.vPeak * this_elFac;
            gain = csum .* this.neScale .* this.vPeak .* this_elFac; % jcl(+)
            
            %Put back into same size as input CTM
            gain = reshape(gain,nRow,nCol);
        end
        
        function setSteer(this,theta,phi,idx_elements)
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
                %applyPhaseQuantization(this)
            end
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
            [azg,elg]=meshgrid(azpts,elpts);
            res = getGain(this,azg,elg,mode);
%             tic;
%             res = zeros(length(elpts),length(azpts),this.precisionClass);
%             hw = waitbar(0,sprintf('Computing Far-Field Gain: %d%% Complete',0));
%             for kaz = 1:length(azpts)
%                 for kel = 1:length(elpts)
%                     res(kel,kaz) = getGain(this,azpts(kaz),elpts(kel),mode);
%                 end
%                 waitbar(kaz/length(azpts),hw,sprintf('Computing Far-Field Gain: %d%% Complete',floor(100*kaz/length(azpts))));
%             end
%             close(hw);
%             tout = toc;
%             fprintf('Total Time: %f sec\nTime per Az/El: %f sec\n',tout,tout/(length(azpts)*length(elpts)));
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
        
        
        function compTaylorWgt(this,NBAR,SLL)
            % Syntax: compTaylorWgt(this,NBAR,SLL)
            % this requires access to MATLAB taylorwin function
            % SLL is input as a negative
            % compute a linear taylor of length = # elements
            taylor_ref_x = linspace(-1,1,this.nElements);
            taylor_ref_y = taylorwin(this.nElements,NBAR,SLL);
            
            % normalize taylor output
            %taylor_ref_y = taylor_ref_y ./ max(abs(taylor_ref_y(:)));
            
            
            % compute normalized range (0 = middle, +/- 1 = edge)
            nrng = this.rng ./ max(abs(this.rng(:)));
            
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