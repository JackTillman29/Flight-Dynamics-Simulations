classdef RcsTool
    %%
    % RcsTool provides a mechanism for viewing, modifying, and exporting
    % ESAMS-formatted signature files.
    %
    % ESAMS RCS signature file format
    %==================================
    %     0  deg  AZ    TAIL
    %    90  deg  AZ    RIGHT WING
    %   180  deg  AZ    NOSE
    %   270  deg  AZ    LEFT WING
    %     0  deg  EL    BELLY
    %    90  deg  EL    WATERLINE
    %   180  deg  EL    TOP
    
    
    %% Properties
    properties
        
        filename = []; % filename being analyzed by the tool
        
        nAz = 0; % number of azimuth points in file
        nEl = 0; % number of elevation points in file
        
        AzDeg = []; % azimuth angle array [deg]
        ElDeg = []; % elevation angle array [deg]
        SqMeter = []; % radar cross section matrix [square meters]
        AzDegMesh = []; % azimuth angle mesh [deg]
        ElDegMesh = []; % elevation angle mesh [deg]
        
    end
    
    properties (Access = private)
        debug_mode = 0; % set this HERE if you want debug code to execute
        fid = []; % file pointer for the RCS file
    end
    
    %% Methods
    % The standard methods for user interaction
    methods
        
        % constructor method
        function this = RcsTool(filename,varargin)
            
            % useful for internal and external temporary class creation
            if(exist('filename'))
                this.filename = filename;
                
                % open file as read-only
                this.fid = fopen(filename,'r');
                
                % load file
                this = this.LoadEsamsSignature();
            else
                this.AzDeg = [-180 0 180];
                this.ElDeg = [0 90 180];
                this.SqMeter = ones(3,3);
                [this.AzDegMesh,this.ElDegMesh]=meshgrid(this.AzDeg,this.ElDeg);
                this.filename = 'notional_0dBsm.rcs';
            end
        end
        
        function this = LoadEsamsSignature(this)
            % Loads in an ESAMS-formatted signature file.
            % 
            % handles blank lines, comments "(", ESAMS-specific keywords (TGTRCS,
            % EXPRCS), files containing lines with multiple values, and
            % files containing lines with single values.
            % 
            % STILL NEED TO FIGURE OUT A WAY IF THERE ARE MULTIPLE AZIMUTH
            % POINTS ON THE SAME LINE OF TEXT.
            
            % ESAMS RCS signature file format
            %     0  deg  AZ    TAIL
            %    90  deg  AZ    RIGHT WING
            %   180  deg  AZ    NOSE
            %   270  deg  AZ    LEFT WING
            %     0  deg  EL    BELLY
            %    90  deg  EL    WATERLINE
            %   180  deg  EL    TOP
            
            load_step = 1;
            
            iaz = 0;
            iel = 0;
            idata = 0;
            new_az_row = 1;
            
            while( ~feof(this.fid) )
                thisline = fgetl(this.fid);
                
                % skip if blank line
                if( isempty(strtrim(thisline)) )
                    continue
                end
                
                % ignore comment lines (begin with "(" character)
                if( strcmp(thisline(1),'(') )
                    continue
                end
                
                % if line is not a comment, but has XXXRCS, then ignore.
                %   TGTRCS, EXPRCS, etc
                if( strcmpi(strtrim(thisline),'rcs') )
                    continue
                end
                
                %==========================================================
                % BELOW THIS LINE, "thisline" SHOULD ONLY BE NUMERIC DATA
                %==========================================================
                
                % get first data
                switch load_step
                    % read the number of Az/El points
                    case 1
                        dat = str2num(thisline);
                        if(length(dat) == 2)
                            this.nAz = dat(1);
                            this.nEl = dat(2);
                            load_step = 2;
                            this.SqMeter = zeros(this.nAz,this.nEl);
                        elseif(length(dat) == 1)
                            if(this.nAz == 0 && dat(1) ~= 0)
                                this.nAz = dat(1);
                            else
                                this.nEl = dat(1);
                                load_step = 2;
                            end
                        end
                    % read the elevation angles
                    case 2
                        tmp = str2num(thisline);
                        iel = iel + 1;
                        
                        % get number of elevation points on this line
                        N = length(tmp);
                        this.ElDeg(iel:(iel + N - 1)) = tmp;
                        iel = iel + N - 1;
                        if(iel == this.nEl)
                            load_step = 3;
                        end
                        
                    % read each azimuth angle and data until end of file
                    case 3
                        tmp = str2num(thisline);
                        % get number of data points on this line
                        N = length(tmp);
                        % if data from previous row is finished being read
                        % in, then "idata" should be zero, and we should
                        % get the next azimuth point
                        if(new_az_row == 1)
                            iaz = iaz + 1;
                            new_az_row = 0;
                            this.AzDeg(iaz) = tmp(1);
                            % subtract one since the first element is the
                            % azimuth angle
                            N = N - 1;
                            % if data is still on this line, then remove
                            % the first element in "tmp", which is the
                            % azimuth angle.
                            if(N > 0)
                                tmp(1) = [];
                            else
                                % check if N == 0, if so, cycle (because no
                                % more data is on this line to be read in)
                                continue
                            end
                        end
                        
                        idata = idata + 1;
                        
                        this.SqMeter(iaz,idata:(idata+N-1)) = tmp;
                        idata = idata + N - 1;
                        
                        % if all data on this row is read in, then reset
                        % "idata" to zero so data in next row can be read
                        % in
                        if(idata == this.nEl)
                            idata = 0;
                            new_az_row = 1;
                        end
                end
            end
            
            % once file is loaded, then ensure it is in:
            %   -180 to 180 deg az format
            %   and
            %   -90 to 90 deg el format
            
            % check for first azimuth of -180
            if ( this.AzDeg(1) ~= -180 )
                if ( (this.AzDeg(1) == 0) & (this.AzDeg(end) == 360)  )
                    warning('This is 0-360 az data. Converting to the -180 to 180 format');
                    j = find(this.AzDeg > 180,1,'first');
                    this.SqMeter = circshift(this.SqMeter,j-1,1);
%                     this.SqMeter = [ ...
%                         this.SqMeter((j-1):end,:); ...
%                         this.SqMeter(2:(j-1),  :) ...
%                         ];
                    this.AzDeg = [ ...
                        this.AzDeg((j-1):end)-360  this.AzDeg(2:(j-1)) ];

                elseif ( (this.AzDeg(1) == 0) & (this.AzDeg(end) == 180) )
                    warning('This is 0-180 az data. Converting to the -180 to 180 format');
                    % need to create mirror image data
                    this.AzDeg = [-this.AzDeg(end:-1:2) this.AzDeg];
                    this.SqMeter  = [this.SqMeter(end:-1:2,:);this.SqMeter];
                else
                    error(['Unknown format. First az = ' num2str(this.AzDeg(1)) ', last az = ' num2str(this.AzDeg(end))]);
                end
            end
            
            % check for first elevation of -90
            %   ** -90 to 90 deg el is not the ESAMS convention; the plotting
            %   in 3d depends on the elevation angles being between 0 and
            %   180 deg (if -90 to 90 is ever wanted in as the convention,
            %   then all that needs to be done in plot3d() is to transform
            %   the ElDeg and ElDegMesh variables from [-90,90] to [0,180]
            %   internally)
%             if ( this.ElDeg(1) ~= -90 )
%                 if ( (this.ElDeg(1) == 0) & (this.ElDeg(end) == 180)  )
%                     warning('This is 0-180 el data. Converting to the -90 to 90 format');
%                     this.ElDeg = this.ElDeg - 90;
% %                     j = find(this.ElDeg > 90,1,'first');
% %                     this.ElDeg = [ ...
% %                         this.ElDeg((j-1):end)-180 this.ElDeg(2:(j-1)) ];
%                 else
%                     error(['Unknown format. First az = ' num2str(this.ElDeg(1)) ', last az = ' num2str(this.ElDeg(end))]);
%                 end
%             end
            
            % produce data in matlab interp2 format
            [this.AzDegMesh,this.ElDegMesh]=meshgrid(this.AzDeg,this.ElDeg);
            
        end
        
        
        function plot3d(this,cutAngleDeg,plotBounds_dB)
            for ircs = 1:length(this)
                if ( nargin < 3 )
                    cutAngleDeg = 0;
                end

                if ( nargin < 2 )
                    a = min(this(ircs).SqMeter(:));
                    b = max(this(ircs).SqMeter(:));
                    plotBounds_dB = 10*log10([a b]);
                    plotBounds_dB(1) = max(plotBounds_dB(1),-60);
                    plotBounds_dB(2) = min(plotBounds_dB(2), 60);
                end
                
                % offset dbsm to have minimum at zero;
                dBsm = 10*log10(this(ircs).SqMeter); %dbsm data (absolute)
                dBsmOffset = min(dBsm(:));
                if(exist('cutPlotBounds_dB'))
                    dBsmOffset = plotBounds_dB(1);
                end
                dBsm_from_zero = dBsm - dBsmOffset;
                dBsm_from_zero = max(dBsm_from_zero,0);
                
                % 3d PointPlot ( in ESAMS RCS frame, x = down, y = right, z = back )
                XGR = dBsm_from_zero.' .* cos(this(ircs).ElDegMesh .* pi/180);
                YGR = dBsm_from_zero.' .* sin(this(ircs).ElDegMesh .* pi/180) .* sin(this(ircs).AzDegMesh .* pi/180);
                ZGR = dBsm_from_zero.' .* sin(this(ircs).ElDegMesh .* pi/180) .* cos(this(ircs).AzDegMesh .* pi/180);

                % 3d PointPlot (in Aircraft Frame, x = Forward, y = Right, z = Down (up for plotting))
                XGA = -ZGR;
                YGA =  YGR;
                ZGA = -XGR;

                % waterline overlay
                ElevationIndex = find(this(ircs).ElDeg >= cutAngleDeg(1),1,'first');
                DBCut = 10*log10(this(ircs).SqMeter(:,ElevationIndex))-dBsmOffset;
                DBCut = reshape(DBCut,1,[]); % put into row vector
                XGRCUT = DBCut .* cos(cutAngleDeg(1) .* pi/180);
                YGRCUT = DBCut .* sin(cutAngleDeg(1) .* pi/180) .* sin(this(ircs).AzDeg .* pi/180);
                ZGRCUT = DBCut .* sin(cutAngleDeg(1) .* pi/180) .* cos(this(ircs).AzDeg .* pi/180);
                XGACUT = -ZGRCUT;
                YGACUT =  YGRCUT;
                ZGACUT =  -XGRCUT;


                %%%figure;
                C = sqrt( ...
                    10.^(0.1*(XGA+dBsmOffset)).^2 + ...
                    10.^(0.1*(YGA+dBsmOffset)).^2 +  ...
                    10.^(0.1*(ZGA+dBsmOffset)).^2);
                C = 10.0 * log10(C);
                C = 0.0*C;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure;
                surf(XGA,YGA,ZGA,C,'EdgeAlpha',0.1,'FaceAlpha',0.4);
                Axes = gca;
                hold on;
                hwl = plot3(XGACUT,YGACUT,ZGACUT,'LineWidth',2);
                hold off;
                set(gca,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1]);
                set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','')
                % set(gca,'Color','k');
                % set(gca,'XTickLabel',abs(get(gca,'XTick'))+dBsmOffset);
                % set(gca,'YTickLabel',abs(get(gca,'YTick'))+dBsmOffset);
                % set(gca,'ZTickLabel',abs(get(gca,'ZTick'))+dBsmOffset);
                view(90,-90);

                if(this.debug_mode)
                    tailstr  = [ ];
                    nosestr  = ['Nose ' ];
                    rwingstr = ['Right Wing ' ];
                    lwingstr = ['Left Wing ' ];
                    topstr   = ['Top ' ];
                    btmstr   = ['Bottom ' ];
                else
                    tailstr  = 'Tail';
                    nosestr  = 'Nose';
                    rwingstr = 'Right Wing';
                    lwingstr = 'Left Wing';
                    topstr   = 'Top';
                    btmstr   = 'Bottom';
                end
                
                hText = text(min(XGA(:)),0,0,        tailstr,'HorizontalAlignment','center');
                hText = [ hText text(max(XGA(:)),0,0,nosestr,'HorizontalAlignment','center') ];

                hText = [ hText text(0,min(YGA(:)),0,rwingstr    ,'HorizontalAlignment','center') ];
                hText = [ hText text(0,max(YGA(:)),0,lwingstr   ,'HorizontalAlignment','center') ];
                hText = [ hText text(0,0,max(ZGA(:)),topstr         ,'HorizontalAlignment','center') ];
                Text = [ hText text(0,0,min(ZGA(:)), btmstr       ,'HorizontalAlignment','center') ];
                view(-90,90);

                %hColorBar = colorbar;
                axis tight;
                %set(hColorBar,'YTickLabel',round(get(hColorBar,'YTick')+dBsmOffset));
                title({ ...
                    'ESAMS RCS Verification Plot', ...
                    this(ircs).filename
                    },'Interpreter','none');
            end
        end
        
        function plotAzCut(this, cutAngleDeg, cutPlotBounds_dB)
            
            compute_plotbounds = 0;
            if(~exist('cutPlotBounds_dB'))
                compute_plotbounds = 1;
            end
            
            for ircs = 1:length(this)    
                figure;
                legendCounter = 1;
                for q = 1:length(cutAngleDeg)
                    [pathstr, namestr] = fileparts(this(ircs).filename);
                    if ( legendCounter == 1 )
                        legCell = {[namestr '@' num2str(cutAngleDeg(q))]};
                    else
                        legCell = {legCell{:},[namestr '@' num2str(cutAngleDeg(q))]};
                    end

                    ElevationIndex = find(this(ircs).ElDeg >= cutAngleDeg(q),1,'first');
                    SqMeterCut = this(ircs).SqMeter(:,ElevationIndex);
                    SqMeterCut = reshape(SqMeterCut,1,[]);
                    % limit so that lines do not go through the origin of the polar plot
                    
                    if(compute_plotbounds)
                        roundNearest5dB = @(x) sign(x).*ceil(abs(x)./5)*5;
                        minmax_cut = 10*log10([min(SqMeterCut) max(SqMeterCut)]);
                        cutPlotBounds_dB = roundNearest5dB(minmax_cut);
                    end
                    
                    SqMeterCut = max(SqMeterCut,10.^(0.1*cutPlotBounds_dB(1)));
                    % prepare a higher-resolution version of non-uniformly
                    % spaced angle samples
                    min_daz = min(abs(diff(this(ircs).AzDeg)));
                    interp_az = linspace(-180,180,ceil(360/min_daz));
                    SqMeterCut = interp1(this(ircs).AzDeg, SqMeterCut, interp_az);
                    
                    if ( nargin >= 2 )
%                         mmpolar_mod(pi/180*(this(ircs).AzDeg+180),10*log10(SqMeterCut),cutPlotBounds_dB);
                        try
                            % polarplot() was added in MATLAB 2016a
                            polarplot(pi/180*(interp_az+180),10*log10(SqMeterCut))
                            set(gca,'RLim',cutPlotBounds_dB)
                            set(gca,'ThetaDir','clockwise')
                            set(gca,'ThetaZeroLocation','top')
                            set(gca,'ThetaLim',[-180 180])
                        catch
                            
                        end

                %           a = [10*log10(SqMeterCut)];

                %           median_0_30_is =  median(a(330:360))
                    else
                %         mmpolar_mod(pi/180*(this(ircs).AzDeg+180),10*log10(SqMeterCut));
                        mmpolar_mod(pi/180*(interp_az+180),10*log10(SqMeterCut),cutPlotBounds_dB);        
%                         pi/180*(this(ircs).AzDeg+180);
%                         length(10*log10(SqMeterCut));
                    end
                    hold on;
                    title({['Azimuth Cut Plot @ ' num2str(cutAngleDeg) ' Elevation'],' '},'Color','k');
                %    hAxes = [hAxes gca];
                    legendCounter=legendCounter+1;
                end
                hold off;

                axisChildren = get(gca,'Children');
                axisColors   = get(gca,'ColorOrder');
                nLines = floor(length(axisChildren)/2);
%                 %DRS: fails here; commenting out for now 
%                 length(this)
%                 try
%                 for k = 1:length(X)
%                 %    set(axisChildren(26+k),'Color',axisColors(k,:));
%                 %      set(axisChildren(k+nLines),'Color',axisColors(k,:));
%                      set(axisChildren(k),'Color',axisColors(k,:));     
%                 end
%                 catch
%                     disp('dan fix please');
%                 end

                legend(legCell,'Interpreter','none','Location','SouthOutside');
            end
        end
        
        
        function plotElCut(this, cutAngleDeg, cutPlotBounds_dB)
            
            compute_plotbounds = 0;
            if(~exist('cutPlotBounds_dB'))
                compute_plotbounds = 1;
            end
            
            for ircs = 1:length(this)    
                figure;
                legendCounter = 1;
                for q = 1:length(cutAngleDeg)
                    [pathstr, namestr] = fileparts(this(ircs).filename);
                    if ( legendCounter == 1 )
                        legCell = {[namestr '@' num2str(cutAngleDeg(q))]};
                    else
                        legCell = {legCell{:},[namestr '@' num2str(cutAngleDeg(q))]};
                    end

                    AzimuthIndex = find(this(ircs).AzDeg >= cutAngleDeg(q),1,'first');
                    SqMeterCut = this(ircs).SqMeter(AzimuthIndex,:);
                    SqMeterCut = reshape(SqMeterCut,1,[]);
                    
                    % limit so that lines do not go through the origin of the polar plot
                    
                    if(compute_plotbounds)
                        roundNearest5dB = @(x) sign(x).*ceil(abs(x)./5)*5;
                        minmax_cut = 10*log10([min(SqMeterCut) max(SqMeterCut)]);
                        cutPlotBounds_dB = roundNearest5dB(minmax_cut);
                    end
                    
                    SqMeterCut = max(SqMeterCut,10.^(0.1*cutPlotBounds_dB(1)));
                    % prepare a higher-resolution version of non-uniformly
                    % spaced angle samples
                    min_del = min(abs(diff(this(ircs).ElDeg)));
                    interp_el = linspace(0,180,ceil(180/min_del));
                    SqMeterCut = interp1(this(ircs).ElDeg, SqMeterCut, interp_el);
                    SqMeterCut = fliplr(SqMeterCut);
                    
                    if ( nargin >= 2 )
                        try
                            % polarplot() was added in MATLAB 2016a
                            polarplot(pi/180*interp_el, 10*log10(SqMeterCut))
                            set(gca,'RLim',cutPlotBounds_dB)
                            set(gca,'ThetaDir','clockwise')
                            set(gca,'ThetaZeroLocation','top')
                            set(gca,'ThetaLim',[-180 180])
                        catch
                            mmpolar_mod(pi/180*interp_el,10*log10(SqMeterCut),cutPlotBounds_dB);
                        end

                %           a = [10*log10(SqMeterCut)];

                %           median_0_30_is =  median(a(330:360))
                    else
                %         mmpolar_mod(pi/180*(this(ircs).AzDeg+180),10*log10(SqMeterCut));
                        mmpolar_mod(pi/180*interp_el,10*log10(SqMeterCut),cutPlotBounds_dB);        
%                         pi/180*(this(ircs).AzDeg+180);
%                         length(10*log10(SqMeterCut));
                    end
                    hold on;
                    title({['Elevation Cut Plot @ ' num2str(cutAngleDeg) ' Azimuth'],' '},'Color','k');
                %    hAxes = [hAxes gca];
                    legendCounter=legendCounter+1;
                end
                hold off;

                axisChildren = get(gca,'Children');
                axisColors   = get(gca,'ColorOrder');
                nLines = floor(length(axisChildren)/2);
%                 %DRS: fails here; commenting out for now 
%                 length(this)
%                 try
%                 for k = 1:length(X)
%                 %    set(axisChildren(26+k),'Color',axisColors(k,:));
%                 %      set(axisChildren(k+nLines),'Color',axisColors(k,:));
%                      set(axisChildren(k),'Color',axisColors(k,:));     
%                 end
%                 catch
%                     disp('dan fix please');
%                 end

                legend(legCell,'Interpreter','none','Location','SouthOutside');
            end
        end
        
        
        function plotSuite(this,cutAngleDeg,plotBounds_dB)
            
            input_cutAngleDeg = 1;
            input_plotBounds_dB = 1;
            if ~exist('cutAngleDeg')
                cutAngleDeg       = 90;
                input_cutAngleDeg = 0;
            end
            if ~exist('plotBounds_dB')
                plotBounds_dB       = [-60 30];
                input_plotBounds_dB = 0;
            end
            
            for ircs = 1:length(this)
                if     (input_cutAngleDeg == 0) & (input_plotBounds_dB == 0)
                    this(ircs).plot;
                    this(ircs).plot3d;
                    this(ircs).plotAzCut(cutAngleDeg,plotBounds_dB)
                elseif (input_cutAngleDeg == 1) & (input_plotBounds_dB == 0)
                    this(ircs).plot;
                    this(ircs).plot3d(cutAngleDeg);
                    this(ircs).plotAzCut(cutAngleDeg,plotBounds_dB)
                elseif (input_cutAngleDeg == 1) & (input_plotBounds_dB == 1)
                    this(ircs).plot;
                    this(ircs).plot3d(cutAngleDeg,plotBounds_dB);
                    this(ircs).plotAzCut(cutAngleDeg,plotBounds_dB)
                end
            end
            
        end
        
        
        function plot(this,varargin)
            % plots a 2D heatmap of the RCS data in file
            %
            % plot(this)
            % plot(this,interp_daz,interp_del)
            % 
            %
            if(nargin == 3)
                interp_daz = varargin{1};
                interp_del = varargin{2};
            else
                interp_daz = 1;
                interp_del = 1;
            end
            
            for ircs = 1:length(this)
                
%                 min_daz = min(diff(this.AzDeg));
%                 min_del = min(diff(this.ElDeg));

%                 % will plot in these arrays
%                 az = -180:interp_daz:180;
%                 el =  -90:interp_del:90;
                
                az = [min(this.AzDegMesh(:)):interp_daz:max(this.AzDegMesh(:))];
                el = [min(this.ElDegMesh(:)):interp_del:max(this.ElDegMesh(:))];

                % ESAMS elevation angles are from 0 to 180, but want to
                % plot this from -90 to 90, therefore add 90 for
                % interpolation points so it matches internal data format
                % of 0 to 180.
                if(min(el) < -45)
                    % if el is < -45, likely it will go from -90 to +90
                    el_offset = 90;
                else
                    el_offset = 0;
                end
                [AZ,EL] = meshgrid(az,el+el_offset);
                
                rcs = interp2(this.AzDegMesh,this.ElDegMesh, this.SqMeter.', AZ,EL);
                
                hfig(ircs) = figure;
    %             hp = pcolor(AZ, EL, 10*log10(RCS))
    %             set(hp,'EdgeColor','none')
                
                % imagesc can be used here BECAUSE the RCS is RESAMPLED on
                % a UNIFORM GRID.
                %    Why imagesc?
                %       Because the data lies in the CENTER of the pixel as
                %       opposed to surf/pcolor which puts the data on the
                %       corners, but colors the face the correct color.
                him(ircs) = imagesc(az, el , 10*log10(rcs));
                hax(ircs) = gca;

                set(hax(ircs),'ydir','normal');
                colorbar
                set(hax(ircs),'XTick',[-180:30:180],'YTick',[-90:15:90])

            end
        end
        
        function globe(this,N,scale)
            minval = -30;
            maxval = 10;
            spanval = maxval - minval;
            if(~exist('scale'))
                pinchval = 0.2;
            else
                pinchval = scale;
            end
            zeroval = minval + 0.5 * spanval;


            [X,Y,Z] = sphere(N);

            % convert X,Y,Z to RCS Az,El
            [AZ,EL,R]=cart2sph(-X,Y,Z);
            EL = EL + pi/2;

            AZ = AZ * 180 / pi;
            EL = EL * 180 / pi;

            GLOBESIG = 10*log10(interp2(this.AzDegMesh,this.ElDegMesh,this.SqMeter.',AZ,EL));
            GLOBESIGLIM = max(GLOBESIG   ,minval);
            GLOBESIGLIM = min(GLOBESIGLIM,maxval);
            SIGDEV = GLOBESIGLIM - zeroval;
            RDEV = SIGDEV * 2 * pinchval / spanval;
            R = R + RDEV;
            [X2,Y2,Z2]=sph2cart(pi/180*(AZ-90),(pi/180*EL-pi/2),R);
            X2 = -X2;

            %figure;
            %hs=surf(X,Y,Z,GLOBESIG);
            %set(gca,'DataAspectRatio',[1 1 1]);
            %set(hs,'EdgeColor','none');

            figure;
            hs=surf(X2,Y2,Z2,GLOBESIG);
            set(gca,'DataAspectRatio',[1 1 1]);
            set(hs,'EdgeColor','none');


            hl  = lightangle(180,0);
            hl2 = lightangle(  0,0);
            material shiny;
            lighting gouraud;
            axis tight;
            set(gca,'XColor','w','YColor','w','ZColor','w','Color','k');
            set(gcf,'Color','k');
            view(-150,30);
            grid off;set(gca,'Visible','off');
            ht_scale = 1.5;
            ht1 = text( ht_scale,0,0,'Right','Color','w');
            ht2 = text(-ht_scale,0,0,'Left' ,'Color','w');
            ht3 = text( 0,ht_scale,0,'Front' ,'Color','w');
            ht4 = text( 0,-ht_scale,0,'Rear' ,'Color','w');
            ht5 = text( 0,0,ht_scale,'Top' ,'Color','w');
            ht6 = text( 0,0,-ht_scale,'Bottom' ,'Color','w');
            set(gca,'DataAspectRatio',[1 1 1]);

            hcb = colorbar; set(hcb,'Color','w')

        end
        
        function [new_rcs]=cpol(hrcs,vrcs)
            % guiding principle: (H + V) / 2.0, unless |H(dB)-V(dB)| > 5dB,
            % then use maximum value
            hdB = 10*log10(hrcs.SqMeter);
            vdB = 10*log10(vrcs.SqMeter);
            diffMatrix = abs(hdB-vdB);
            useMaxMatrix = diffMatrix > 5.0;
            MaxMatrix = max(hrcs.SqMeter,vrcs.SqMeter);
            new_rcs = hrcs; % make a copy from hrcs
            new_rcs.SqMeter = ((hrcs.SqMeter + vrcs.SqMeter) .* 0.5 .* ~useMaxMatrix) + ...
                               useMaxMatrix .* MaxMatrix;
        end
        
        
        function [varargout] = bugsplat_atmos(S,R0,R0E,Pt,DCY,Gt,Lt,FreqGHz,xgrid,ygrid,altitude,varargin)
            % bugsplat(S,R0,xgrid_m,ygrid_m,altitude_ft)
            % bugsplat(S,R0,-50000:1000:50000,-50000:1000:50000,alt_ft)
            % if(nargin == 2 && isstruct(R0)) % structure input
            %     input = R0; % map 2nd input to input structure
            %     Pt = input.Pt;
            %     DCY = input.DCY;
            %     R0 = input.R0;
            %     R0E = elevation (deg) of R0 number
            %     Gt = input.Gt;
            %     Lt = input.Lt;
            %     FreqGHz = input.FreqGHz;
            %     xgrid = input.xgrid;
            %     ygrid = input.ygrid;
            %     altitude = input.altitude;
            %     structureFormat = 1;
            % elseif(nargin == 2)
            %     error('Structure argument issue.');
            % end
            % 
            %if(nargin ~= 1)
                % FreqGHz is used for atmospheric loss calculations
                load AtmosphericTable.mat;
                % dB per kilometer
                % lookUpExample = interp1(InputGHz,OutputdBPerKm,10.0)

                % compute Ro constants
                Pt = 10*log10(Pt*DCY);

                % atmospheric loss coefficient dB per meter
                alpha = interp1(InputGHz,OutputdBPerKm,FreqGHz) * 1.0e-3; % dBperKm to dBperm
            %end
                outputPointValues = [];
            %if(~structureFormat)
                % setup default values
                input.jam_bandwidth_hz = [];
                input.threat_bandwidth_hz = [];
                input.jammer_range = [];  % empty indicates on-board
                input.jammer_rel_pos = []; % relative downrange,crossrange
                input.jammer_threat_gain_dbi = [];
                input.jammer_coherence_loss_db = [];
                input.jammer_max_power = [];
                input.rcs_plot_bounds_db = [];
                input.output_units = 'km';
                input.output_sf = 1.0e-3;
                input.cancel_attenDb = 0;
                input.maxDBCOLOR = 999;
                input.minDBCOLOR = -999;
                input.J_to_S_SolverDb = 0;
                input.pitchAngleDeg = 0;
                input.rcsVsRange = [];
                input.debug = 0;
                if(exist('Gt'))
                    input.Gr    = Gt;
                else
                    Gr = [];
                end
                input.pulsewidths = [];
                input.specific_pulse = [];
                input.pulse_comp_deltas = [];
                input.noAtmos = 0;

                input.outputPoints = [];
                if(nargin == 1)
                    input.R0 = [];
                    input.Pt = [];
                    input.DCY = [];
                    input.Gt = [];
                    input.Lt = [];
                    input.FreqGHz = [];
                    input.xgrid = [];
                    input.ygrid = [];
                    input.altitude = [];
                    varargout{1} = input;
                    return
                end



                % input processing
                %disp(length(varargin))
                for k = 1:2:length(varargin)
                    thisString = varargin{k};
                    thisValue  = varargin{k+1};
                    %disp(thisString);
                    switch thisString
                        case 'JamBandwidthHz'
                            input.jam_bandwidth_hz = thisValue;
                        case 'ThreatBandwidthHz'
                            input.threat_bandwidth_hz = thisValue;
                        case 'JammerRangeMeters'
                            input.jammer_range = thisValue;
                        case 'JammerRelPos'
                            input.jammer_rel_pos = thisValue;
                        case 'JammerThreatGainDbi'
                            input.jammer_threat_gain_dbi = thisValue;
                        case 'JammerCoherenceLossDb'
                            input.jammer_coherence_loss_db = thisValue;
                        case 'JammerMaxPowerWatts'
                            input.jammer_max_power = 10*log10(thisValue);
                        case 'JammerCancellerLossDb'
                            input.cancel_attenDb = thisValue;
                        case 'RcsPlotLimitsDbsm'
                            input.rcs_plot_bounds_db = thisValue;
                        case 'MaxDbWDisplay'
                            input.maxDBCOLOR = thisValue;
                        case 'MinDbWDisplay'
                            input.minDBCOLOR = thisValue;
                        case 'J_to_S_SolverDb'
                            input.J_to_S_SolverDb = thisValue;
                        case 'PitchAngleDeg'
                            input.pitchAngleDeg = thisValue;
                        case 'RcsVsRange'
                            input.rcsVsRange = thisValue;
                        case 'debug'
                            input.debug = thisValue;
                        case 'Gr'
                            input.Gr = thisValue;
                        case 'noAtmos'
                            input.noAtmos = thisValue;
                        case 'OutputPoints'
                            input.outputPoints = thisValue;
                        case 'Pulsewidths'
                            input.pulsewidths = thisValue;
                        case 'SpecificPulse'
                            input.specific_pulse = thisValue;
                        case 'PulseCompDeltas'
                            input.pulse_comp_deltas = thisValue;
                        case 'OutputUnits'
                            testval = thisValue;
                            switch testval
                                case 'm'
                                    input.output_units = 'm';
                                    input.output_sf = 1.0;
                                case 'km'
                                    input.output_units = 'km';
                                    input.output_sf = 1.0e-3;
                                case 'nmi'
                                    input.output_units = 'nmi';
                                    input.output_sf = 1.0 / 1852;
                                otherwise
                                    error(['Unknown OutputUnits option: ' testval]);
                            end
                        otherwise
                            error(['Unknown option: ' num2str(thisString) ]);
                    end
                end

                % input checking
                if(~isempty(input.jammer_range) & isempty(input.jammer_threat_gain_dbi))
                    error('You must specify a threat gain (dBi) in the direction of the jammer if you specify a jammer range');
                end
                if(~isempty(input.jam_bandwidth_hz) & isempty(input.threat_bandwidth_hz))
                    error('You must specify a threat bandwidth if you specify a jammer bandwidth');
                end
                if(isempty(input.jam_bandwidth_hz) & ~isempty(input.threat_bandwidth_hz))
                    error('You must specify a jammer bandwidth if you specify a threat bandwidth');
                end

                if(isempty(input.jammer_threat_gain_dbi))
                    input.jammer_threat_gain_dbi = input.Gr;
                end

                if(input.noAtmos)
                    alpha = 0;
                end

                if(input.debug)
                    disp([ 'Atmospheric Loss: ' num2str(1e3*alpha) ' dB/km']);
                end

                if(Pt < 20.0)
                    warning(['Radar power input is ' num2str(Pt) 'dBW. Seems low. User input value should be in Watts.']);
                end

            %end


            [X,Y,Z]=meshgrid(xgrid,ygrid,0.3048*altitude);

            R = sqrt(X.^2 + Y.^2 + Z.^2);
            RdB = 10.0*log10(R);  % ranges of all grid points from site (dBmeters)

            %atmLossDb = alpha .* (R - R0);
            % 
            % X2 = X(1:10:end,1:10:end);
            % Y2 = Y(1:10:end,1:10:end);
            % Z2 = Z(1:10:end,1:10:end);

            %**************************************************
            %***************ATMLOS RETURNS ATTENUATION IN LINEAR (ONE-WAY)
            %**************************************************

            % for k1 = 1 : size(X2,1)
            %     for k2 = 1 : size(X2,2)
            %         atmLossDb_temp(k1,k2) = atmlos_byp(S, ...
            %             X2(k1,k2), ...
            %             Y2(k1,k2), ...
            %             Z2(k1,k2), ...
            %             .3 / FreqGHz);
            %     end
            % end
            %atmLossDbTgt = interp2(X2,Y2,atmLossDb_temp,X,Y);
            atmLossDbTgt = this.atmlos_byp( ...
                        X, ...
                        Y, ...
                        Z, ...
                        .3 / FreqGHz);
            % Compute reference atmospheric attenuation at the r0
            atmLossDb_R0 = this.atmlos_byp( ...
                sqrt(R0^2 - (0.3048*altitude)^2), ...
                0, ...
                0.3048*altitude, ...
                .3 / FreqGHz, 1, R0E);


            % Compute grid-based atmospheric loss relative to R0 anchor point
            atmLossDb = (10*log10(1.0./atmLossDbTgt) - 10*log10(1.0./atmLossDb_R0));

            if(input.noAtmos == 1)
                atmLossDb = 0;
                disp('Zeroing out atmospheric losses due to input request.');
            end


            if(~isempty(input.pulsewidths))
                % Generate Pulsewidth matrix based on range
                PW_IDX = 0*R;
                for k = length(input.pulsewidths) : -1 : 1
                    thisBlankRange = input.pulsewidths(k) * 150.0;
                    idxs = find(R > thisBlankRange);
                    PW_IDX(idxs) = k;
                end
                if(~isempty(input.specific_pulse))
                    PW_IDX(:) = input.specific_pulse;
                end

                % Generate pulse compression delta matrix
                % **** Assumption is that largest (1st) pulsewidth is associated with
                %      input R0 term.
                R0_Delta = 0*R;
                JamCohLoss = 0*R;

                for k = length(input.pulsewidths) : -1 : 1
                    thisPulsewidth = k;
                    idxs = find(PW_IDX==k);
                    R0_Delta(idxs) = input.pulse_comp_deltas(k);
                end

                % Generate pulse compression loss jammer matrix
                if(isempty(input.jammer_coherence_loss_db))

                elseif(length(input.jammer_coherence_loss_db) == 1)
            %        disp(['Assuming same coherence loss for all waveforms!']);
            %        for k = length(input.pulsewidths) : -1 : 1
            %             thisPulsewidth = k;
            %             idxs = find(PW_IDX==k);
            %             JamCohLoss(idxs) = input.jammer_coherence_loss_db; % - input.pulse_comp_deltas(k);
            %        end 
                    error('this does not work yet. spec all losses per waveform.');
                else
                    for k = length(input.pulsewidths) : -1 : 1
                        thisPulsewidth = k;
                        idxs = find(PW_IDX==k);
                        JamCohLoss(idxs) = input.jammer_coherence_loss_db(k) - input.pulse_comp_deltas(k);
                    end 
                end



                if(input.debug)
                    figure;
                    htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), PW_IDX' );
                    xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                    colormap(jet(length(input.pulsewidths)));
                    caxis([1 length(input.pulsewidths)+1]);
                    hc=colorbar;

                    title('Pulsewidth Selection');
                    axis square;grid on;
                    set(gca,'YDir','normal');
                    set(hc,'YTick',1.5:1:(length(input.pulsewidths)+0.5),'YTickLabel', ...
                        input.pulsewidths);

                    figure;
                    htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), R0_Delta' );
                    xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                    %colormap(jet(length(input.pulsewidths)));
                    %caxis([1 length(input.pulsewidths)+1]);
                    hc=colorbar;

                    title('Pulse Compression Delta');
                    axis square;grid on;
                    set(gca,'YDir','normal');
                    %set(hc,'YTick',1.5:1:(length(input.pulsewidths)+0.5),'YTickLabel', ...
                    %    input.pulsewidths);

                    figure;
                    htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), JamCohLoss' );
                    xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                    %colormap(jet(length(input.pulsewidths)));
                    %caxis([1 length(input.pulsewidths)+1]);
                    hc=colorbar;

                    title('Jammer Compression Loss');
                    axis square;grid on;
                    set(gca,'YDir','normal');
                end 
            else
                if(~isempty(input.jammer_coherence_loss_db))
                    JamCohLoss = 0*R+input.jammer_coherence_loss_db;
                else
                    JamCohLoss = 0*R;
                end
            end

            % figure;
            % imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), atmLossDb' );
            % colorbar;
            % axis square;grid on;
            % set(gca,'YDir','normal');


            % Generate Jammer Range Array
            if(~isempty(input.jammer_range)) % fixed jammer range
                tmprcs = RcsTool();
                Rjam = 10*log10(0.0*R + input.jammer_range);
                H = 0.3048 * altitude;
                jamAtmLoss = tmprcs.atmlos_byp(sqrt(input.jammer_range^2 - H^2),0,H, .3 / FreqGHz);
                atmLossDbJam = 0*R + 10*log10(jamAtmLoss);
            elseif(isempty(input.jammer_rel_pos)) % jammer range = self protect
                Rjam = 10*log10(R);
                atmLossDbJam = atmLossDbTgt;
            else
                XJAM = X-(sign(X).*input.jammer_rel_pos(1));
                YJAM = Y-(sign(Y).*input.jammer_rel_pos(2));
                ZJAM = Z;
                Rjam = 10.0*log10(sqrt(XJAM.^2 + YJAM.^2 + ZJAM.^2));

            %     XJAM2 = XJAM(1:10:end,1:10:end);
            %     YJAM2 = YJAM(1:10:end,1:10:end);
            %     ZJAM2 = ZJAM(1:10:end,1:10:end);
            %     
            %     for k1 = 1 : size(XJAM2,1)
            %         for k2 = 1 : size(XJAM2,2)
            %             atmLossDbJam_temp(k1,k2) = atmlos_byp( ...
            %                 XJAM2(k1,k2), ...
            %                 YJAM2(k1,k2), ...
            %                 ZJAM2(k1,k2), ...
            %                 .3 / FreqGHz);
            %         end
            %     end
            %     
            %     atmLossDbJam = interp2(XJAM2,YJAM2,atmLossDbJam_temp,XJAM,YJAM);

                atmLossDb = this.atmlos_byp( ...
                            XJAM, ...
                            YJAM, ...
                            ZJAM, ...
                            .3 / FreqGHz);

            end
            RjamLin = 10.^(0.1*Rjam);

            %====================================


            %atmLossDbJam = alpha .* (RjamLin - R0);
            %disp('fix this...');

            if(input.debug)
                figure;
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), atmLossDb' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Atmospheric Loss Applied to Target Location (One-Way)');
                axis square;grid on;
                set(gca,'YDir','normal');

                figure;
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), atmLossDbJam' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Atmospheric Loss Applied to Jammer (rel R0 range)');
                axis square;grid on;
                set(gca,'YDir','normal');

                figure;
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), input.output_sf*RjamLin' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Jammer Linear Range');
                axis square;grid on;
                set(gca,'YDir','normal');

            end

            %bugsplat_atmos(sig,7*1852,7*1852,240,1,26,0,16.8,-50000:500:50000,-50000:500:50000,0,'OutputUnits','nmi')

            % compute azimuth and elevation grid values
            % These angles define the radar LOS to the target
            ZPRIME = Z.*sind(input.pitchAngleDeg) + Z.*cosd(input.pitchAngleDeg);
            XPRIME = X.*cosd(input.pitchAngleDeg) - X.*sind(input.pitchAngleDeg);
            EL = 90-180/pi*asin(ZPRIME./max(R,1));

            EL = 0*EL + 90 + 10; % 90 is waterline

            AZ = -180/pi*atan2(Y,-XPRIME);

            if(input.debug)
                figure;
                subplot(1,2,1);
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), AZ' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Azimuth RCS Lookup');
                axis square;grid on;
                set(gca,'YDir','normal');
                subplot(1,2,2);
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), EL' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Elevation RCS Lookup');
                axis square;grid on;
                set(gca,'YDir','normal');
            end


            %EL = 90-180/pi*asin(Z./R);  % add 90 to get to ESAMS reference (180=top, 0 = bottom)
            %AZ = -180/pi*atan2(Y,-X);

            % Lookup the RCS for each grid point based upon the actual RCS deck and the
            % AZ,EL angles
            if(isempty(input.rcsVsRange))
                sigma_grid = 10*log10(griddata(S.ElMeshDeg,S.AzMeshDeg,S.SqMeter',EL,AZ));
            else
                sigma_grid = 10*log10(interp1(input.rcsVsRange(:,1),input.rcsVsRange(:,2),R));
            end

            sigma_plus_atm_grid = sigma_grid - 2.0 * atmLossDb;

            if(~isempty(input.pulsewidths))
                sigma_plus_atm_grid = sigma_plus_atm_grid + R0_Delta;
            end

            if(input.debug)
                figure;
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), sigma_plus_atm_grid' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Equivalent RCS dBsm (wv changes and L_atm applied two way at this location)');
                axis square;grid on;
                set(gca,'YDir','normal');
            end



            % The detection range along any LOS will be exactly the 0dBsm detection
            % range degraded by the presented RCS. If the RCS goes down by 3dB, the
            % detection range goes down by 3dB/4 or 0.75 dBmeters
            DetectionRangeDb = 10*log10(R0)+0.25*sigma_plus_atm_grid;
            DetectionRangeDeficitDb = DetectionRangeDb - RdB;

            % If the detection range deficit is positive, this means that the target
            % can be detected at that grid point. The power needed is related to the
            % range by P = 20*log10(R), or PdB = RdB * 2.0
            %PowerNeededSelfProtect = DetectionRangeDeficitDb .* 2.0;

            % Mismatch and Coherence Loss
            SignalLoss = 0*R; %input.jammer_coherence_loss_db;
            if(~isempty(input.jam_bandwidth_hz))
                SignalLoss = SignalLoss + 10*log10(input.jam_bandwidth_hz/input.threat_bandwidth_hz);
            end
            % Derate the self protect power as appropriate
            PowerNeeded = JamCohLoss + (Pt + Gt + input.Gr + sigma_plus_atm_grid + 2*Rjam - Lt) - (input.jammer_threat_gain_dbi + 10*log10(4*pi) + 4*RdB) + SignalLoss + input.cancel_attenDb + atmLossDbJam + input.J_to_S_SolverDb;
            %                      *  *  *  *  *   *    *    *  *   *
            % Verified: PJ(REQ) = PT+GT+GR-LT+RCS+2RJ-4PIDB-4R-GRJ+LJ+2LA-LAJ
            %
            %  Note: RCS here = (sigma - 2*atmLossDb) , atmLossDb is positive on the
            %  denominator, thus (sigma + 2*atmLossDb).
            %        LAJ here is positive on the denominator, thus LAJ = -LAJ

            if(isempty(input.jammer_max_power))
                input.jammer_max_power = max(PowerNeeded(:));
            end

            if(input.debug)
                figure;
                subplot(1,2,1);
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), SignalLoss' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title('Jammer Waveform & Spectral Signal Loss');
                axis square;grid on;
                set(gca,'YDir','normal');

                subplot(1,2,2);
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), PowerNeeded' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colorbar;
                title(['Jammer Power Required, dBW, for ' num2str(input.J_to_S_SolverDb) ' J/S']);
                axis square;grid on;
                set(gca,'YDir','normal');
            end



            % build the alpha mask for max detection
            alpha_mask_det = ~(DetectionRangeDeficitDb < 0);

            % =======  Presented RCS Plot ========== %
            figure;
            hImage = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), sigma_grid' );
            if(nargout > 0)
                varargout{1} = sigma_grid';
                varargout{2} = alpha_mask_det';
            end

            set(hImage,'AlphaData',alpha_mask_det');
            disp('skipping alphadata');
            xlabel(['Downrange (' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
            if(~isempty(input.rcs_plot_bounds_db))
                caxis(input.rcs_plot_bounds_db);
            end
            colorbar;
            title({'\fontsize{16}Presented RCS (dBsm)',['R_{0}: ' num2str(input.output_sf*R0) ' ' input.output_units ', Alt: ' num2str(altitude) 'ft']});
            axis square;grid on;
            set(gca,'YDir','normal');


            % Check to see if a blank range mask is needed
            BlankMask = 0.0 * R + 1;
            if(~isempty(input.pulsewidths))
                if(~isempty(input.specific_pulse))
                    [Rmin,Rmax,RminP,RmaxP] = UsefulRangeLfm(1e-6*input.pulsewidths(input.specific_pulse),10e-3,0.7);
                    %BlankMask = (R > Rmin);
                    %set(hImage,'AlphaData',alpha_mask_det' .* BlankMask');
                end
            end




            % ====== Required Power Plot =========== %
            if(~isempty(input.jammer_max_power))
                % Mask: Only data where power needed is less than max power and the
                % platform was detectable
                %alpha_mask = ~(PowerNeeded > input.jammer_max_power | DetectionRangeDeficitDb < 0);
                alpha_mask_power_sat = (PowerNeeded > input.jammer_max_power);


                alpha_mask = ((.8*alpha_mask_power_sat + .2*alpha_mask_det)).*alpha_mask_det;
                iTooHigh = find(PowerNeeded > input.maxDBCOLOR);
                alpha_mask(iTooHigh) = 0.0;

                figure;
                %hImage2 = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), ~alpha_mask2');set(hImage2,'AlphaData',0.1*alpha_mask2');
                %hold on;
                hImage = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), PowerNeeded' );
                hold off;
                set(hImage,'AlphaData',alpha_mask');
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                if(input.maxDBCOLOR == 999)
                    %caxis([min(PowerNeeded(:)) max(PowerNeeded(:))]);
                    maxc = max(PowerNeeded(:));
                else
                    %caxis([min(PowerNeeded(:)) input.maxDBCOLOR]);
                    maxc = input.maxDBCOLOR;
                end

                if(input.minDBCOLOR == -999)
                    minc = min(PowerNeeded(:));
                else
                    minc = input.minDBCOLOR;
                end
                caxis([minc maxc]);

                colorbar;
                %title({'Predicted Detection Contour',['R_{0}: ' num2str(R0) ', Alt: ' num2str(altitude./.3048) 'ft'],'Color Coding = Presented RCS'});
                axis square;grid on;
                set(gca,'YDir','normal');

                % == title
                if(~isempty(input.jammer_range))
                    jloc = ['Jamming Offboard @ Range of ' num2str(input.output_sf*input.jammer_range) ' ' input.output_units];
                else
                    jloc = 'Jamming Onboard';
                end
                title({'Power Required to Cover Platform',jloc,...
                    ['Threat Bandwidth (kHz): ' num2str(1e-3*input.threat_bandwidth_hz) ', Jam Bandwidth (kHz): ' num2str(1e-3*input.jam_bandwidth_hz)] ...
                    ['Threat Gain to Jammer (dB): ' num2str(input.jammer_threat_gain_dbi)], ...
                    ['Canceller Loss (dB): ' num2str(input.cancel_attenDb)], ...
                    ['Jam Waveform Coherence Loss (dB): ' num2str(input.jammer_coherence_loss_db)]});

                if(exist('Rmin'))
                    [xpc,ypc]=circle_patch(0,0,input.output_sf.*Rmin,100);
                    %hpblank = patch(xpc,ypc,'r');
                    %set(hpblank,'EdgeColor','none','FaceColor',0.25*[1 1 1],'FaceAlpha',1.0);
                end
            end

            if( length(input.pulsewidths) > 1 )
                figure;
                htemp = imagesc( input.output_sf*unique(Y), input.output_sf*unique(X), PW_IDX' );
                xlabel(['Downrange(' input.output_units ')']);ylabel(['Crossrange(' input.output_units ')']);
                colormap(jet(length(input.pulsewidths)));
                caxis([1 length(input.pulsewidths)+1]);
                hc=colorbar;
                set(htemp,'AlphaData',alpha_mask_det');

                title('Pulsewidth Selection');
                axis square;grid on;
                set(gca,'YDir','normal');
                set(hc,'YTick',1.5:1:(length(input.pulsewidths)+0.5),'YTickLabel', ...
                    input.pulsewidths);
            end

            if(~isempty(input.outputPoints))
                optr = 1;
                for k = 1 : 2 : length(input.outputPoints)
                    thisx = input.outputPoints(k);
                    thisy = input.outputPoints(k+1);
                    this_sigma = interp2(xgrid,ygrid,sigma_grid,thisx,thisy);
                    this_det   = interp2(xgrid,ygrid,alpha_mask_det,thisx,thisy);
                    this_power = interp2(xgrid,ygrid,PowerNeeded,thisx,thisy);

                    if(this_det == 0)
                        this_power = nan;
                    end    
                    fprintf('%8.3f,%8.3f: Sigma = %8.3f dBsm, Power = %8.3f dBW\n',input.output_sf*thisx,input.output_sf*thisy,this_sigma,this_power);
                    outputPointValues(optr,:) = [this_sigma this_power];
                    optr = optr + 1;
                end
                varargout{3} = outputPointValues;
            end
        end
        
    end
    
    
    methods
        function [aloss,dB_per_km] = atmlos_byp(this,X,Y,Z,WAVEL,bypass,bypass_elevation_deg)

            %% Constants and Code Setup

            SPDLGT = 2.99792458e8;
            REARTH = 8497000;
            WAVEL = SPDLGT / 9.8e9;
            xspan   = 300000;
            yspan   = 200000;
            xstep   = 20000;
            alt_agl = 9144;
            ant_agl = 5;

            % [X,Y,Z] = meshgrid( ...
            %     -xspan:xstep:xspan, ...
            %     -yspan:xstep:yspan, ...
            %     alt_agl-ant_agl);

            FREQS = [.1e9,.2e9,.3e9,.6e9,1.e9,3.e9,10.e9];
            ELEVS = [0.,.5,1.,2.,5.,10.];

            AA = [.2739, .1881, .1605, .1031, .07371, .04119;
                .6848, .5533, .4282, .3191, .2158,  .1017;
                1.199, .9917, .7498, .5186, .3029,  .1522;
                2.210, 1.830, 1.314, .9499, .4724,  .2512;
                2.758, 2.177, 1.798, 1.168, .5732,  .3007;
                3.484, 2.592, 1.964, 1.345, .6478,  .3408;
                4.935, 3.450, 2.601, 1.718, .9130,  .4420]';

            BB = [   .008648, .008644, .01106, .01723, .02313, .04076;
                .008648, .008644, .01104, .01374, .02213, .04886;
                .006837, .008795, .01110, .01474, .03116, .05360;
                .008499, .009737, .01221, .01623, .03677, .07204;
                .01030,  .01223,  .01163, .01831, .03927, .08056;
                .009745, .01225,  .01455, .02055, .04500, .08280;
                .00999,  .01340,  .01620, .02240, .03750, .08470 ]';



            %% Round Earth computations
            R  = sqrt(X.^2 + Y.^2 + Z.^2);
            GR = sqrt(X.^2 + Y.^2       );

            if(exist('bypass','var'))
                ELEV=ones(size(X,1),size(X,2)) .* bypass_elevation_deg;
            else
                E1 = (REARTH^2+R.^2 - (REARTH+Z).^2) ./ (2*REARTH.*R);
                ELEV = 180/pi * acos(E1) - 90.0;
                ELEV(GR == 0) = 90.0;
            end

            RNM = min(R ./ 1852,300); % 300 nmi limit
            FREQ = min(max(SPDLGT / WAVEL,FREQS(1)),FREQS(end));
            ELEV = min(max(ELEV,ELEVS(1)),ELEVS(end));
            FREQ_ARRAY = 0*ELEV + FREQ;

            %% Interpolation Section (all input data is bounded)
            % I,J indexing maps to elevation,frequency

            A = interp2(FREQS,ELEVS,AA,FREQ(1:end),ELEV(1:end));
            B = interp2(FREQS,ELEVS,BB,FREQ(1:end),ELEV(1:end));
            A = reshape(A,size(X,1),size(X,2));
            B = reshape(B,size(X,1),size(X,2));
            LOSSDB = A.*(1-exp(-B.*RNM));%two way log
            LOSS = 10.^(0.1*LOSSDB);%two way lin
            ALOSSDB = 5*log10(LOSS);%one way log
            aloss=10.^(-0.1*ALOSSDB);%one way lin
            dB_per_km = ALOSSDB ./ R;

        end

        
        
    end
    
    methods
        function writeEsamsSignatureFile(this,filename,descriptionString)
        % Syntax:
        %   writeEsamsSignatureFile(RCSData,filename,descriptionString)    
            maxPerLine = 5;
            fid = fopen(filename,'w');
            if(iscell(descriptionString))
                for istr = 1:length(descriptionString)
                    fprintf(fid,'%s\n',['( *** ' descriptionString{istr}]);
                end
            else
                fprintf(fid,'%s\n',['( *** ' descriptionString]);
            end
            fprintf(fid,'%s\n','( This file was automatically generated from WriteEsamsSignatureFile.m');
            fprintf(fid,'TGTRCS %d %d',[length(this.AzDeg) length(this.ElDeg)]);
            fprintf(fid,'%s','    ');
            fprintf(fid,'\n%12.2f%12.2f%12.2f%12.2f%12.2f',this.ElDeg);
            for k = 1:length(this.AzDeg)
                fprintf(fid,'%s\n','');
                fprintf(fid,'%8.2f',this.AzDeg(k));
                fprintf(fid,'\n%12.3e%12.3e%12.3e%12.3e',this.SqMeter(k,:));
            end
            fclose(fid);
        end
        
        function rcsdat = getRcsLookup(this,azDeg,elDeg,interp_method)
            % Returns the interpolated RCS values at user-input az/el points
            % (default method is 'nearest')
            % 
            % The user can input scalars, vectors, and matrices for the
            % az/el points.
            if ~exist('interp_method')
                interp_method = 'nearest';
            end
            
            rcsdat = interp2(this.AzDegMesh, this.ElDegMesh, this.SqMeter, ...
                azDeg, elDeg, interp_method);
            
        end
        
    end

end