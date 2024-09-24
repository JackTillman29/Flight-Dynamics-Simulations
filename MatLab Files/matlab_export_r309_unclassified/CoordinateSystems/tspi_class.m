% Assumptions: 
% Input units will be decimal deg lat, lon
% Altitude (m)
% Attitude (YPR) can be either deg or rad (class field).
% Output is always in meters & radians.
classdef tspi_class < handle
    properties
        filename;
        hasYPR = 0;
        tableData;
        nPlats = 0;
        nRuns = 0;
        platIDs = {};
        modelMode = 1;  % 1-esams 3dof, 2-sup 3dof, 3-esams 6dof
        
        units_roll  = 'deg';
        units_pitch = 'deg';
        units_yaw   = 'deg';
        
        % These are the default search strings. You may override them.
        searchString_wgs84_lat = 'LAT84';
        searchString_wgs84_lon = 'LON84';
        searchString_wgs84_alt = 'ALT_MSL';
        searchString_time      = 'TSEC';
        searchString_roll      = 'ROLL';
        searchString_pitch     = 'PT';
        searchString_yaw       = 'HD';
        searchString_run       = 'RUN';
    end
    methods
        function this = tspi_class(infile,mode)
            this.filename = infile;
            this.tableData = readtable(infile);
            this.reprocess_file();
            % manage consecutive rows with same time
            this.cleanup_table;
        end
        
        function cleanup_table(this)
            time    = getField(this,this.searchString_time);
            keepIdx = find(diff(time)>0)+1;
            this.tableData = this.tableData([1; keepIdx],:);
        end
        
        function reprocess_file(this)
            % check for roll
            oc = this.findOccurences(this.searchString_roll);
            if(oc > 0)
                this.hasYPR = 1;
            end
            
            % check # of LAT84s for # of platforms
            [oc,oc_names] = this.findOccurences(this.searchString_wgs84_lat);
            this.nPlats = length(oc);
            this.platIDs = this.replaceString(oc_names,this.searchString_wgs84_lat,'');
            if(~isempty( this.getField(this.searchString_run) ) )
                this.nRuns = length(unique(getField(this,this.searchString_run)));
            else
                this.nRuns = 1;
            end
            if(exist('mode','var'))
                this.modelMode = mode;
            end
        end
        
        function showHeader(this)
            disp(this.tableData(1,:));
        end
        
        function [oc,ocn,ocnk] = findOccurences(this,str)
            % returns occurence idx, and cell array of occurence names
            oc = [];
            ocn = {};
            for k = 1 : length(this.tableData.Properties.VariableNames)
                test_str = this.tableData.Properties.VariableNames{k};
                if(strfind(test_str,str))
                    oc = [oc k];
                    ocn = {ocn{:} test_str};
                end
            end
        end
        
        function res = replaceString(this,in,sub_from,sub_to)
            res = {};
            for k = 1 : length(in)
                res = {res{:}, strrep(in{k},sub_from,sub_to)};
            end
        end
        
        function LLA84 = getLLA(this,varargin)
            %   getLLA(this, platform, run_number)
            %   getLLA(this, varargin):
            %       'String-Value' pairs:
            %          'platform'  platform_number
            %          'run'       run_number
            for k = 1:2:length(varargin)
                switch lower(varargin{k})
                    case 'platform'
                        platform = varargin{k+1};
                    case 'run'
                        run_num = varargin{k+1};
                end
            end
            
            if(~ischar(platform))
                platform = num2str(platform);
            end
            
            idxLat = findOccurences(this,[this.searchString_wgs84_lat platform]);
            if(isempty(idxLat))
                error(['Could not locate header: ' this.searchString_wgs84_lat platform]);
            end
            
            idxLon = findOccurences(this,[this.searchString_wgs84_lon platform]);
            if(isempty(idxLon))
                error(['Could not locate header: ' this.searchString_wgs84_lon platform]);
            end
            
            idxAlt = findOccurences(this,[this.searchString_wgs84_alt platform]);
            if(isempty(idxAlt))
                error(['Could not locate header: ' this.searchString_wgs84_alt platform]);
            end
            
            LLA84 = table2array(this.tableData(:,[idxLat idxLon idxAlt]));
            if(exist('run_num','var'))
                if(this.nRuns > 1)
                    [~,runidx]=this.getField(this.searchString_run);
                    run_num_rows = find(table2array(this.tableData(:,runidx)) == run_num);
                elseif( (this.nRuns == 1) &&  (run_num == 1))
                    run_num_rows = 1:size(this.tableData,1);
                else
                    error(['Could not locate run # ' num2str(run_num) ' in file!']);
                end
                LLA84 = LLA84(run_num_rows,:);
            elseif(this.nRuns > 1)
                warning([this.filename ' has multiple runs!! You are not filtering by run number!']);
            end
            
        end
        
        
        function setLLA(this,LLA84,varargin)
            %   setLLA(this, varargin):
            %       'String-Value' pairs:
            %          'platform'  platform_number
            %          'run'       run_number
            for k = 1:2:length(varargin)
                switch lower(varargin{k})
                    case 'platform'
                        platform = varargin{k+1};
                    case 'run'
                        run_num = varargin{k+1};
                end
            end
            
            if(~ischar(platform))
                platform = num2str(platform);
            end
            
            idxLat = findOccurences(this,[this.searchString_wgs84_lat platform]);
            if(isempty(idxLat))
                error(['Could not locate header: ' this.searchString_wgs84_lat platform]);
            end
            
            idxLon = findOccurences(this,[this.searchString_wgs84_lon platform]);
            if(isempty(idxLon))
                error(['Could not locate header: ' this.searchString_wgs84_lon platform]);
            end
            
            idxAlt = findOccurences(this,[this.searchString_wgs84_alt platform]);
            if(isempty(idxAlt))
                error(['Could not locate header: ' this.searchString_wgs84_alt platform]);
            end
            
            if(exist('run_num','var'))
                if(this.nRuns > 1)
                    [~,runidx]=this.getField(this.searchString_run);
                    run_num_rows = find(table2array(this.tableData(:,runidx)) == run_num);
                elseif( (this.nRuns == 1) &&  (run_num == 1))
                    run_num_rows = 1:size(this.tableData,1);
                else
                    error(['Could not locate run # ' num2str(run_num) ' in file!']);
                end
                this.tableData(run_num_rows,idxLat) = table(LLA84(:,1));
                this.tableData(run_num_rows,idxLon) = table(LLA84(:,2));
                this.tableData(run_num_rows,idxAlt) = table(LLA84(:,3));
            elseif(this.nRuns > 1)
                warning([this.filename ' has multiple runs!! You are not filtering by run number!']);
                % check to see that input data is same length as the
                % overall table:
                if size(LLA84,1) == size(this.tableData,1)
                    this.tableData(:,idxLat) = table(LLA84(:,1));
                    this.tableData(:,idxLon) = table(LLA84(:,2));
                    this.tableData(:,idxAlt) = table(LLA84(:,3));
                end
            end
            
            
            
        end
        
        
        
        function [dat,datidx] = getField(this,field_name,run_num)
            % pass in EXACT field name to extract
            [idxdat,idxstr] = findOccurences(this,field_name);
            dat = [];
            datidx = [];
            for k = 1 : length(idxdat)
                if(strcmp(idxstr{k},field_name))
                    dat = table2array(this.tableData(:,idxdat(k)));
                    if(exist('run_num','var'))
                        [~,runidx]=this.getField(this.searchString_run);
                        if(~isempty(runidx))
                            run_num_rows = find(table2array(this.tableData(:,runidx)) == run_num);
                            dat = dat(run_num_rows,:);
                        end
                    end
                    datidx = idxdat(k);
                end
            end
        end
        
        
        
        function YPR = getYPR(this,varargin)
            %   getYPR(this, varargin):
            %       'String-Value' pairs:
            %          'platform'  platform_number
            %          'run'       run_number
            
            % check if the first input is a string and not a number
            for k = 1:2:length(varargin)
                switch lower(varargin{k})
                    case 'platform'
                        platform = varargin{k+1};
                    case 'run'
                        run_num = varargin{k+1};
                    otherwise
                        error('unknown input string')
                end
            end
            
            if(~ischar(platform))
                platform = num2str(platform);
            end
            if(this.hasYPR)
                idxRoll = findOccurences(this,[this.searchString_roll platform]);
                if(isempty(idxRoll))
                    error(['Could not locate header: ' this.searchString_roll platform]);
                end

                idxPitch = findOccurences(this,[this.searchString_pitch platform]);
                if(isempty(idxPitch))
                    error(['Could not locate header: ' this.searchString_pitch platform]);
                end

                idxYaw = findOccurences(this,[this.searchString_yaw platform]);
                if(isempty(idxYaw))
                    error(['Could not locate header: ' this.searchString_yaw platform]);
                end

                YPR = table2array(this.tableData(:,[idxYaw idxPitch idxRoll]));
                if(exist('run_num','var'))
                    if(this.nRuns > 1)
                        [~,runidx]=this.getField(this.searchString_run);
                        run_num_rows = find(table2array(this.tableData(:,runidx)) == run_num);
                    elseif( (this.nRuns == 1) &&  (run_num == 1))
                        run_num_rows = 1:size(this.tableData,1);
                    else
                        error(['Could not locate run # ' num2str(run_num) ' in file!']);
                    end
                    YPR = YPR(run_num_rows,:);
                elseif(this.nRuns > 1)
                    warning([this.filename ' has multiple runs!! You are not filtering by run number!']);
                end
            else
                error('Detected no attitude information in this file.');
            end
            
            % Unit modification
            switch lower(this.units_yaw)
                case 'rad'
                    YPR(:,1) = 1.0 * YPR(:,1);
                case 'deg'
                    YPR(:,1) = pi/180 * YPR(:,1);
                otherwise
                    error(['Incorrect yaw angle units: ' this.units_yaw]);
            end
            
            switch lower(this.units_pitch)
                case 'rad'
                    YPR(:,2) = 1.0 * YPR(:,2);
                case 'deg'
                    YPR(:,2) = pi/180 * YPR(:,2);
                otherwise
                    error(['Incorrect pitch angle units: ' this.units_pitch]);
            end
            
             switch lower(this.units_roll)
                case 'rad'
                    YPR(:,3) = 1.0 * YPR(:,3);
                case 'deg'
                    YPR(:,3) = pi/180 * YPR(:,3);
                otherwise
                    error(['Incorrect roll angle units: ' this.units_roll]);
            end
            
        end
        
        function setYPR(this,YPR,varargin)
            %   setYPR(this, varargin):
            %       'String-Value' pairs:
            %          'platform'  platform_number
            %          'run'       run_number
            
            % check if the first input is a string and not a number
            for k = 1:2:length(varargin)
                switch lower(varargin{k})
                    case 'platform'
                        platform = varargin{k+1};
                    case 'run'
                        run_num = varargin{k+1};
                    otherwise
                        error('unknown input string')
                end
            end
            
            if(~ischar(platform))
                platform = num2str(platform);
            end
            
            % Unit modification
            switch lower(this.units_yaw)
                case 'rad'
                    YPR(:,1) = 1.0 * YPR(:,1);
                case 'deg'
                    YPR(:,1) = pi/180 * YPR(:,1);
                otherwise
                    error(['Incorrect yaw angle units: ' this.units_yaw]);
            end
            
            switch lower(this.units_pitch)
                case 'rad'
                    YPR(:,2) = 1.0 * YPR(:,2);
                case 'deg'
                    YPR(:,2) = pi/180 * YPR(:,2);
                otherwise
                    error(['Incorrect pitch angle units: ' this.units_pitch]);
            end
            
             switch lower(this.units_roll)
                case 'rad'
                    YPR(:,3) = 1.0 * YPR(:,3);
                case 'deg'
                    YPR(:,3) = pi/180 * YPR(:,3);
                otherwise
                    error(['Incorrect roll angle units: ' this.units_roll]);
             end
            
            % apply data to table
            if(this.hasYPR)
                idxRoll = findOccurences(this,[this.searchString_roll platform]);
                if(isempty(idxRoll))
                    error(['Could not locate header: ' this.searchString_roll platform]);
                end

                idxPitch = findOccurences(this,[this.searchString_pitch platform]);
                if(isempty(idxPitch))
                    error(['Could not locate header: ' this.searchString_pitch platform]);
                end

                idxYaw = findOccurences(this,[this.searchString_yaw platform]);
                if(isempty(idxYaw))
                    error(['Could not locate header: ' this.searchString_yaw platform]);
                end

                if(exist('run_num','var'))
                    if(this.nRuns > 1)
                        [~,runidx]=this.getField(this.searchString_run);
                        run_num_rows = find(table2array(this.tableData(:,runidx)) == run_num);
                    elseif( (this.nRuns == 1) &&  (run_num == 1))
                        run_num_rows = 1:size(this.tableData,1);
                    else
                        error(['Could not locate run # ' num2str(run_num) ' in file!']);
                    end
%                     this.tableData(run_num_rows,[idxYaw idxPitch idxRoll]) = YPR;
                    this.tableData(run_num_rows,idxYaw)   = table(YPR(:,1));
                    this.tableData(run_num_rows,idxPitch) = table(YPR(:,2));
                    this.tableData(run_num_rows,idxRoll)  = table(YPR(:,3));
                elseif(this.nRuns > 1)
                    warning([this.filename ' has multiple runs!! You are not filtering by run number!']);
                    if size(YPR,1) == size(this.tableData,1)
%                         this.tableData(:,[idxYaw idxPitch idxRoll]) = YPR;
                        this.tableData(:,idxYaw)   = table(YPR(:,1));
                        this.tableData(:,idxPitch) = table(YPR(:,2));
                        this.tableData(:,idxRoll)  = table(YPR(:,3));
                    end
                end
            else
                error('Detected no attitude information in this file.');
            end
            
        end
        
        function [tVec,outData] = interpData(this,inTime,inData,dt)
            % create time vector
            tVec = inTime(1):dt:inTime(end);
            outData = interp1(inTime,inData,tVec);
            
        end
        
        function writeESAMS_TSPI(this,platID,run_num,wind)
            outdat = runEstimator(this,platID,run_num,wind);
            [~,basename]=fileparts(this.filename);
            fid = fopen([basename '_run' num2str(run_num) '_' platID '.tgfp'],'w');
            fprintf(fid,'%18.6f %18.6f %18.6f %18.6f %18.6f %18.6f %18.6f %18.6f %18.6f %18.6f\n',outdat);
            fclose(fid);
        end
        
        function varargout = dumpESAMS_ECEF(this,varargin)
            %addpath('I:\MATLAB\DigitalFilters');
            
            if(length(varargin) > 0)
                for kk = 1:2:length(varargin)
                    switch lower(varargin{kk})
                        case {'platform','plat'}
                            plat_array = varargin{kk+1};
                        case {'run','run_num','runnum'}
                            run_array = varargin{kk+1};
                    end
                end
            else
                if(isempty(unique(this.getField(this.searchString_run))'))
                    run_array = 1;
                else
                    run_array = unique(this.getField(this.searchString_run))';
                end
                plat_array = 1:this.nPlats;
            end
            
            ifile = 0;
            for isel = 1:length(run_array)
                krun  = run_array(isel);
                kplat = plat_array(isel);
%             % loop over runs in the file
%             for krun = run_array
%                 % loop on number of platforms
%                 for i = plat_array
                    
                    ifile = ifile + 1;
                    
                    % determine starting column of lat/lon/alt data
                    %firstColumn = 5 + (i-1)*3;
                    clear wgs84Data wgs84ypr;
                    % gather WGS-84 data for estimator tool
                    wgs84Data(:,2:4) = this.getLLA('platform',num2str(kplat),'run',krun);%   dataTable.data(:,3);
                    wgs84Data(:,1)   = this.getField(this.searchString_time,krun);
                    wgs84Data(:,1)   = wgs84Data(:,1) - wgs84Data(1,1);
                    tout = wgs84Data(:,1);
                    %wgs84Data(:,3) = dataTable.data(:,firstColumn+1);
                    %wgs84Data(:,4) = dataTable.data(:,firstColumn+2);
                    
                    %[wgsout,tout]=this.despikeLLA(wgs84Data(:,1),wgs84Data(:,2:4));
                    wgsout = wgs84Data(:,2:4);
                    tout   = wgs84Data(:,1);
                    
                    
                    wgs84Data(:,1) = tout;
                    wgs84Data(:,2:4) = wgsout;
                    wgs84ypr = this.getYPR('platform',num2str(kplat),'run',krun);
                    ecefData = convertGeodeticWGS84LatLonAltToECEF(wgs84Data(:,2) , wgs84Data(:,3) , wgs84Data(:,4)*.3048 );
                    
                    ecefFilterData = StateEstimator(tout,ecefData,1e-6,1,1e6,1e6,1e-5);
                    ecefYPR = zeros(length(tout),3);
                    
                    % Now we need ECEF YPR
                    % This function computes the uvw basis on the wgs84 ellipsoid under
                    % the aircraft in ECEF coordinates
                    for ri = 1 : length(wgs84Data(:,1))
                        ECEF2NED_LTP = ComputeWGS84_ECEF_Basis( ...
                            wgsout(ri,1), ...    % WGS84 Lattitude in degrees
                            wgsout(ri,2) );      % WGS84 Longitude in degrees
                        %if(ri == 1)
                        %    assignin('base','ECEF2NED_LTP',ECEF2NED_LTP);
                        %end
                        % compute 3x3 from LTP to body
                        T_LTP2B = RotationMatrixYPR( ...
                            wgs84ypr(ri,1), ...
                            wgs84ypr(ri,2), ...
                            wgs84ypr(ri,3));
                        
                        %if(ri == 1)
                        %    assignin('base','T_LTP2B',T_LTP2B);
                        %end
                        
                        T = T_LTP2B*ECEF2NED_LTP;
                        
                        [y,p,r]=Tmx321_to_YPR(T);
                        ecefYPR(ri,:) = [y p r];
                        
                    end
                    outdat = [tout ecefFilterData.xhat(:,1:6) ecefYPR(:,[2 3 1])];
                    %assignin('base','outdat',outdat);
                    % Write output data
                    [pathTo,name_noext,~]=fileparts(this.filename);
                    output_file{ifile} = [pathTo '\' name_noext '_run' num2str(krun) '_plat' num2str(kplat) '.tgfpe'];
                    if(nargout == 0)
                        % only write this message to command window when
                        % not outputting the filenames as a cell array of
                        % strings.
                        disp(['Writing ' output_file{ifile}]);
                    end
                    fid = fopen(output_file{ifile},'w');
                    %for ii = 1:length(enu_pos_out(:,1))
                    % negate roll angle for ESAMS TGFP output
                    fprintf(fid,'4\n!### ESAMS Earth Centered Earth Fixed (ECEF) Entity Flight Path ###\n');
                    fprintf(fid,'!### t(s),WGS84 Lat (deg), WGS84 Lon (deg), WGS84 Alt (ft), WGS84 Yaw (rad), WGS84 Pitch (rad), WGS84 Roll (rad), x(m),y(m),z(m),vx(m/s),vy(m/s),vz(m/s),p(rad),r(rad),y(rad) ###\n');
                    fprintf(fid,'!### *** NOTE: WGS84 values are for reference only!! They are not used by ESAMS and may be zeros. *** ###\n');
                    fprintf(fid,'!### ***       For ECEF flight paths generated using tspi_class.m, they will be included for convenience ###\n');
                    fprintf(fid,'!### Attitude angles are associated with Euler 321 sequence ###\n');
                    fprintf(fid,'!### from ECEF to A/C body frame, right handed, F-R-D ###\n');
                    fprintf(fid,'!### Generated on: %s %s %s\n',datestr(now),'by',getenv('USERNAME'));
                    fprintf(fid, ...
                        '%12.6f %12.6f %12.4f %8.2f %8.4f %8.4f %8.4f %10.2f %10.2f %10.2f %10.3f %10.3f %10.3f %11.7f %11.7f %11.7f\n',...
                        [ outdat(:,1) wgs84Data(:,1:3) wgs84ypr outdat(:,2:end) ]');
                    
                    %end
                    fclose(fid);
                    
                    simdisYaw = outdat(:,10);
                    simdisRoll = outdat(:,9);
                    simdisPitch = outdat(:,8);
                    
                    outdat2 = outdat;
                    outdat2(:,5:7) = [simdisYaw simdisPitch simdisRoll];
                    outdat2(:,8:10) = outdat(:,5:7);
                    
                    fid = fopen([output_file{ifile} '.asi'],'w');
                    fprintf(fid,'Version     22\nRefLLA      37.0 -115.0 0.0\nCoordSystem "ECEF"\n');
                    fprintf(fid,'PlatformID 1\nPlatformName 1 "Platform 1"\nPlatformIcon 1 "f-15c_eagle"\n');
                    fprintf(fid,'PlatformData 1 %8.2f %10.2f %10.2f %10.2f %10.3f %10.3f %10.3f %11.7f %11.7f %11.7f\n',outdat2');
                    fclose(fid);
                    
%                 end
%             end
            end
            
            if(nargout > 0)
                varargout{1} = output_file;
            end
            
        end
        
    end
end
