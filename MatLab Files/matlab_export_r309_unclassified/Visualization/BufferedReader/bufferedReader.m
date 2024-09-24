classdef bufferedReader < handle
% --------- bufferedReader help ------------
%
% Syntax: bufferedReader(filename,headBytes,footBytes,pageSize,precision)
%
% Supported precision: single,double,ubit14
    properties
        headBytes;
        footBytes;
        recBytes;   %bytes per record?
        recBytesWhdr;
        recLength;  %2, one for I, one for Q
        valBytes;
        precision;
        filename;
        fileBytes; %total bytes in file
        fid;

        %number of 4-byte records (2 for I, 2 for Q) in file
        %nRecords is fileBytes/4
        nRecords;  
        
        remRecords;
        pageStartRec; %(current page-1)*pageSize + 1
        pageStopRec;  %currentpage * pageSize
        pageStartByte;  %
        pageStopByte;
        pageMemory;
        pageSize;
        pageRunning = 0;
        fileloc = -1;  %location of file pointer, measured in bytes
        endian;
        fread_machinefmt;
        headBlockBytes;
        pageOverlap;
        partialShift_numRecords;
    end
    
    properties (Constant=true)
        % shiftPage argument enumeration constants (for easier code reading)
        enum_fast = 1;
        enum_load = 0;
    end
    
    methods
        function this = bufferedReader(filename,headBytes,footBytes,pageSize,precision,varargin)
            this.filename = filename;
            this.headBytes = headBytes;
            this.footBytes = footBytes;
            this.pageSize = pageSize;
            this.precision = precision; % single or double
            this.headBlockBytes = 0;
            this.pageOverlap = 0;
            this.endian = 'native';
            this.fread_machinefmt = 'native';
            if(length(varargin) == 1)
            % little endian:        (memory)
                %       32-bit integer  |   .  |
                %       0A 0B 0C 0D     |   .  |
                %        |  |  |  |     |   .  |
                %        |  |  |  --->  |  0D  |
                %        |  |  ------>  |  0C  |
                %        |  --------->  |  0B  |
                %        ------------>  |  0A  |
                %                       |  .   |
                %
                % big endian:
                %                       (memory)
                %                       |   .  |  32-bit integer
                %                       |   .  |     0A 0B 0C 0D
                %                       |   .  |     |  |  |  |
                %                       |  0A  |  <---  |  |  |
                %                       |  0B  |  <------  |  |
                %                       |  0C  |  <---------  |
                %                       |  0D  |  <------------
                %                       |  .   |
                switch varargin{1}
                    case 'BigEndian'
                        this.endian = 'BigEndian';
                        this.fread_machinefmt = 'ieee-be';
                    case 'LittleEndian'
                        this.endian = 'LittleEndian';
                        this.fread_machinefmt = 'ieee-le';
                    case 'BigEndian64'
                        this.endian = 'BigEndian64';
                        this.fread_machinefmt = 'ieee-be.l64';
                    case 'LittleEndian64'
                        this.endian = 'LittleEndian64';
                        this.fread_machinefmt = 'ieee-le.l64';
                    otherwise
                        this.endian = 'native';
                        this.fread_machinefmt = 'native';
                end
            end
            if(length(varargin) > 0)
                varstep = 2;
                for ivar = 1:varstep:length(varargin)
                    varstep = 2;
                    if(strcmpi(varargin{ivar},'endian'))
                        % see above diagram for explanation of "endianness"
                        switch varargin{ivar+1}
                            case 'BigEndian'
                                this.endian = 'BigEndian';
                                this.fread_machinefmt = 'ieee-be';
                            case 'LittleEndian'
                                this.endian = 'LittleEndian';
                                this.fread_machinefmt = 'ieee-le';
                            case 'BigEndian64'
                                this.endian = 'BigEndian64';
                                this.fread_machinefmt = 'ieee-be.l64';
                            case 'LittleEndian64'
                                this.endian = 'LittleEndian64';
                                this.fread_machinefmt = 'ieee-le.l64';
                            otherwise
                                this.endian = 'native';
                                this.fread_machinefmt = 'native';
                        end
                    elseif(strcmpi(varargin{ivar},'HeaderBlockBytes'))
                        this.headBlockBytes = varargin{ivar+1};
                    elseif(strcmpi(varargin{ivar},'recordLength'))
                        this.recLength = varargin{ivar+1};
                    elseif(strcmpi(varargin{ivar},'pageOverlap'))
                        this.pageOverlap = varargin{ivar+1};
                    elseif(strcmpi(varargin{ivar},'partialShiftNumRecords'))
                        this.partialShift_numRecords = varargin{ivar+1};
                    end
                end
            end
            
            switch(lower(precision))
                case 'single'
                    this.valBytes = 4;
                case 'double'
                    this.valBytes = 8;
                case 'ubit14'
                    this.valBytes = 14/8;
                case 'int14'
                    this.valBytes = 14/8;
                case 'uint16'
                    this.valBytes = 2;
                case 'int16'
                    this.valBytes = 2;
                case 'int8'
                    this.valBytes = 1;
                otherwise
                    errordlg(['Unknown precision passed! ' precision]);
                    return;
            end
            
            % initialize
            
            % check for file
            if(~exist(this.filename))
                errordlg(['Could not locate ' this.filename]);
                error(['Could not locate ' this.filename]);
            end
            disp([' Loading: ' this.filename]);
            
            % open file and assess
            this.fid = fopen(this.filename,'rb');
            
            if(this.fid < 1)
                errordlg(['Error opening file ' this.filename]);
            end
            
            % get file size
            fseek(this.fid,0,'eof');
            this.fileBytes = ftell(this.fid) - this.headBlockBytes;
            disp(['  Total Bytes: ' num2str(this.fileBytes)]);
            fseek(this.fid,0,'bof');
            
            % Read header (sanity check)
            if(this.headBytes ~= 0)
                % uint32 is data type specific to FORTRAN's write function
                % TODO: make this a variable input (possibly enumerate for FORTRAN)
                this.recBytes = fread(this.fid,1,'uint32') - this.headBlockBytes;
                this.recBytesWhdr = this.recBytes + this.headBytes + this.footBytes;
                this.recLength = this.recBytes/this.valBytes;
            else
                if(~isempty(this.recLength))
                    this.recBytes = this.valBytes*this.recLength;% - this.headBlockBytes;
                    this.recBytesWhdr = this.recBytes;
                else
                    % KDS - can't do this!! A binary file may not always
                    % contain the # of bytes as the first value. (e.g. BSAG
                    % data).
                    if(this.headBlockBytes == 0 & this.footBytes == 0)
                        warning('Because no header or footer, assuming record length of 1!!');
                        this.recBytes = this.valBytes;
                        this.recBytesWhdr = this.recBytes;
                        this.recLength = 1;
                    else
                        this.recBytes = fread(this.fid,1,'uint32');
                        this.recBytesWhdr = this.recBytes;
                        this.recLength = this.recBytes/this.valBytes;
                    end
                end
%                 this.recLength = ;
%                 this.recBytes = this.valBytes;% - this.headBlockBytes;
%                 this.recBytesWhdr = this.valBytes;
%                 this.recLength = 1;
            end
            disp(['  Record size (bytes): ' num2str(this.recBytes)]);
            disp(['  Record size (vals) : ' num2str(this.recLength)]);
            fseek(this.fid,this.headBlockBytes+1,'bof');
            this.fileloc = this.headBlockBytes+1;
            
            this.nRecords = double(this.fileBytes / this.recBytesWhdr);
            disp(['  # Records in file : ' num2str(this.nRecords)]);
            
            if(~isempty(this.partialShift_numRecords))
                % if partial page shift was specified in # of samples to
                % shift each page by, then compute the corresponding page
                % overlap that yields this shift
                this.pageOverlap = (this.pageSize - this.partialShift_numRecords)/this.pageSize;
            end
            
            % if page size is not specified "[]", use full file
            if(isempty(pageSize))
                pageSize = this.nRecords;
                this.pageSize = this.nRecords;
            end
            
            % make sure page isn't larger than data
            if(pageSize > this.nRecords)
                pageSize = this.nRecords;
                this.pageSize = pageSize;
                %warndlg('Page size exceeds data available. Truncating!');
                warning('Page size exceeds data available. Truncating!');
            end
            
            % setup page memory
            bitloc = strfind(precision,'14');
            if(bitloc)
                s1 = precision(1);
                if(s1 == 'u')
                    s1 = '';
                end
                
                se = num2str(precision((bitloc+3):end));
                
                prec_store = [ s1 'int' num2str(2^ceil(log2(14))) ];
                    
                
                this.pageMemory = zeros(pageSize,this.recLength,'int16');
                warning(['using ' prec_store ' storage class.']);
            else
                this.pageMemory = zeros(pageSize,this.recLength,lower(precision));
            end
            
            resetPage(this);
            
        end
        
        function resetPage(this)
            % because we're at the beginning
            this.remRecords = this.nRecords;
            
            %this.pageStartRec  = 1;
            %this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
            %this.pageStartByte = this.headBytes;
            %this.pageStopByte  = this.pageStartByte - this.headBytes + this.pageSize*this.recBytesWhdr;
            
            this.pageStartRec  = 1;
            this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
            this.pageStartByte = 0 + this.headBlockBytes + 1; % couldn't get
            %working for MIDAS.. seems like we need it though.
%             this.pageStartByte = 0;

            %Original Line
            %this.pageStopByte  = this.pageStartByte + int64(this.pageSize*this.recBytesWhdr) - 1 + this.headBlockBytes;
            %new line (CTM, 9/20/17)
            this.pageStopByte  = this.pageStartByte + (this.pageSize*this.recBytesWhdr) - 1 + this.headBlockBytes;

            frewind(this.fid);
            fseek(this.fid,this.headBlockBytes,'cof');
            
            %original line
            %this.fileloc = int64(this.headBlockBytes)+1;
            %new line (CTM, 9/20/17)
            this.fileloc = (this.headBlockBytes)+1;
            this.pageRunning = 0;
            nextPage(this);
            this.pageRunning = 1;
        end
        
        function resetPageFast(this)
            % because we're at the beginning
            this.remRecords = this.nRecords;
            
            %this.pageStartRec  = 1;
            %this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
            %this.pageStartByte = this.headBytes;
            %this.pageStopByte  = this.pageStartByte - this.headBytes + this.pageSize*this.recBytesWhdr;
            
            this.pageStartRec  = 1;
            this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
            this.pageStartByte = 0 + this.headBlockBytes + 1; % couldn't get
            %working for MIDAS.. seems like we need it though.
%             this.pageStartByte = 0;

            %Original Line
            %this.pageStopByte  = this.pageStartByte + int64(this.pageSize*this.recBytesWhdr) - 1 + this.headBlockBytes;
            %new line (CTM, 9/20/17)
            this.pageStopByte  = this.pageStartByte + (this.pageSize*this.recBytesWhdr) - 1 + this.headBlockBytes;

            frewind(this.fid);
            fseek(this.fid,this.headBlockBytes,'cof');
            
            %original line
            %this.fileloc = int64(this.headBlockBytes)+1;
            %new line (CTM, 9/20/17)
            this.fileloc = (this.headBlockBytes)+1;
            this.pageRunning = 0;
            shiftPage(this,0,enum_fast);
            this.pageRunning = 1;
        end
        
        
        
        function nextPage(this)
            if(this.pageOverlap == 0)
                shiftPage(this,0,this.enum_load); % 0 indicates "next"
            else
                partialShiftPage(this,0);
            end
        end
        
        
%         function jumpToPage(this,n)
%             disp(ftell(this.fid));
%             if(n <= this.nRecords)
%                 this.resetPage();
%                 for k = 1 : (n-1)
%                     this.nextPage();
%                 end
%             else
%                 error('Out of bounds');
%             end
%             disp(ftell(this.fid));
%         end
        
        function jumpToPage(this,newPage)

            if(newPage <= this.nRecords)
                oldPage = (this.pageStartRec-1)/this.pageSize + 1;
                if round(oldPage)~=oldPage
                    error('oldPage is not an integer');
                end;

                if newPage<oldPage  

                    if newPage==1
                        this.resetPage();
                    else
                        this.resetPageFast(); %go back to beginning
                        for k = 1:newPage-2
                            shiftPage(this,0,this.enum_fast); %going fowards
                        end;
                        nextPage(this);
                    end;

                elseif newPage>oldPage %going forward
                    for k = 1:(newPage-oldPage-1)
                        shiftPage(this,0,this.enum_fast); %going forward
                    end;
                    nextPage(this);
                else
                    return;
                end
            else
                error('Out of bounds');
            end
        end;
        
        function prevPage(this)
            shiftPage(this,1,0); % 1 indicates "previous"
        end
        
        function shiftPage(this,direction,fast)
            fullpage = 1;
            nRecRead = single(this.pageSize);
            
            if(this.pageRunning == 1)
                % data is populated so bounds checking must be performed
                
                
                if(direction == 0) % going forwards
                    if(this.remRecords >= this.pageSize)
                        % all good.
                        this.pageStartRec  = this.pageStopRec + 1;
                        this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
                        this.pageStartByte = (this.pageStartRec - 1) * this.recBytesWhdr+double(this.headBlockBytes);
                        this.pageStopByte  = this.pageStartByte + this.pageSize*this.recBytesWhdr - 1;
                        fullpage = 1;
                    elseif(this.remRecords > 0)
                        % can only show some of the data
                        this.pageStartRec  = this.pageStopRec + 1;
                        this.pageStopRec   = this.pageStartRec + this.remRecords - 1;
                        this.pageStartByte = (this.pageStartRec - 1) * this.recBytesWhdr+double(this.headBlockBytes);
                        this.pageStopByte  = this.pageStartByte + this.remRecords*this.recBytesWhdr-1;
                        fullpage = 0;
                        nRecRead = this.remRecords;
                    else
                        disp('End of file');
                        return;
                    end
                else % going backwards
                    %this.remRecords = max(this.remRecords,0);
                    recordsInFront = this.nRecords - (this.pageSize+this.remRecords);
                    if(recordsInFront >= this.pageSize)
                        % all good.
                        %this.pageStartRec  = this.pageStopRec + 1 - 2*this.pageSize;
                        firstRewind = this.pageStopRec-this.pageStartRec+1;
                        % TODO TODO TODO
                        this.pageStartRec  = this.pageStopRec - (firstRewind+this.pageSize-1);
                        this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
                        this.pageStartByte = (this.pageStartRec - 1) * this.recBytesWhdr;
                        this.pageStopByte  = this.pageStartByte + this.pageSize*this.recBytesWhdr-1;
                        fullpage = 1;
                        
                        %fseek(this.fid,-((this.pageSize+firstRewind)*this.recBytesWhdr-this.footBytes),'cof');
                        fseek(this.fid,this.pageStartByte-this.fileloc,'cof');
                        this.fileloc = ftell(this.fid);
                    else
                        disp('Beginning of file');
                        return;
                    end
                end
            end
            
            % check for alignment - at this point, we assume that the
            % file position is at a record boundary
            if( this.fileloc ~= (this.pageStartByte) )
                error(['Byte misalignment detected in ' this.filename]);
            end
            
            % advance one header
            fseek_result = fseek(this.fid,this.headBytes,'cof');
            this.fileloc = ftell(this.fid);
            
            
            if(fast == this.enum_fast)
                % shift page pointers only w/o bulk read
                fseek(this.fid,2*double(nRecRead*this.recLength),'cof');
            else
                % bulk read in data
                pstring = [num2str(this.recLength) '*' this.precision];
                temp = fread( ...
                    this.fid, ...
                    double(nRecRead*this.recLength), ...
                    pstring, ...
                    this.footBytes+this.headBytes, ...
                    this.fread_machinefmt);
            end
            
            % back up one header (fread skip option skips footer and next header)
            % in order to get back to record boundary, this must be done
            fseek(this.fid,-this.headBytes,'cof');
            this.fileloc = ftell(this.fid);
            
            this.pageMemory(1:nRecRead,:) = reshape(temp, ...
                this.recLength, ...
                nRecRead).';
                        
            this.pageMemory((nRecRead+1):end,:) = 0;
            if(direction == 0)
                this.remRecords = this.remRecords-this.pageSize;
            else
                this.remRecords = this.remRecords+this.pageSize;
            end
        end
        
        function partialShiftPage(this,direction)
            fullpage = 1;
            nRecRead = single(this.pageSize);
            
            jumpBackBytes = 0;
            if(this.pageRunning == 1)
                % data is populated so bounds checking must be performed
                
                
                if(direction == 0) % going forwards
                    if(this.remRecords >= this.pageSize)
                        % all good.
%                         disp(['1:  ' num2str(this.pageStartRec) '    ' num2str(this.pageStopRec)])
%                         disp(['1:  ' num2str(this.pageStartByte) '    ' num2str(this.pageStopByte)])
                        
                        this.pageStartRec  = this.pageStartRec + floor(this.pageSize*(1-this.pageOverlap)*this.recLength)/this.recLength;
                        this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
%                         this.pageStartByte = ceil((this.pageStartRec - 1) * this.recBytesWhdr)+double(this.headBlockBytes);
                        this.pageStartByte = round((this.pageStartRec-1) * this.recBytesWhdr * this.valBytes)+double(this.headBlockBytes);
%                         this.pageStartByte = (this.pageStartRec - 1) * double(this.recBytesWhdr)+double(this.headBlockBytes);

%                         after  = this.pageStartByte;
%                         even = (after-before)/this.recBytes
                        this.pageStopByte  = this.pageStartByte + this.pageSize*this.recBytesWhdr - 1;
%                         disp(['2:  ' num2str(this.pageStartRec) '    ' num2str(this.pageStopRec)])
                        disp(['2:  ' num2str(this.pageStartByte) '    ' num2str(this.pageStopByte) '   delta:' num2str(this.pageStopByte-this.pageStartByte)])
%                         disp(['fileloc:  ',num2str(this.fileloc)])
                        pause(0.1)
%                         this.pageStartRec  = this.pageStopRec + 1;
%                         this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
                        %                         this.pageStopByte  = this.pageStartByte + this.pageSize*this.recBytesWhdr - 1;
                        fullpage = 1;
                    elseif(this.remRecords > 0)
                        % can only show some of the data
                        this.pageStartRec  = this.pageStopRec + 1;
                        this.pageStopRec   = this.pageStartRec + this.remRecords - 1;
                        this.pageStartByte = (this.pageStartRec - 1) * this.recBytesWhdr+double(this.headBlockBytes);
                        this.pageStopByte  = this.pageStartByte + this.remRecords*this.recBytesWhdr-1;
                        fullpage = 0;
                        nRecRead = this.remRecords;
                    else
                        disp('End of file');
                        return;
                    end
                else % going backwards
                    %this.remRecords = max(this.remRecords,0);
                    recordsInFront = this.nRecords - (this.pageSize+this.remRecords);
                    if(recordsInFront >= this.pageSize)
                        % all good.
                        %this.pageStartRec  = this.pageStopRec + 1 - 2*this.pageSize;
                        firstRewind = this.pageStopRec-this.pageStartRec+1;
                        % TODO TODO TODO
                        this.pageStartRec  = this.pageStopRec - (firstRewind+this.pageSize-1);
                        this.pageStopRec   = this.pageStartRec + this.pageSize - 1;
                        this.pageStartByte = (this.pageStartRec - 1) * this.recBytesWhdr;
                        this.pageStopByte  = this.pageStartByte + this.pageSize*this.recBytesWhdr-1;
                        fullpage = 1;
                        
                        %fseek(this.fid,-((this.pageSize+firstRewind)*this.recBytesWhdr-this.footBytes),'cof');
                        fseek(this.fid,this.pageStartByte-this.fileloc,'cof');
                        this.fileloc = ftell(this.fid);
                    else
                        disp('Beginning of file');
                        return;
                    end
                end
            end
            
            % check for alignment - at this point, we assume that the
            % file position is at a record boundary
%             if( this.fileloc ~= (this.pageStartByte) )
%                 error(['Byte misalignment detected in ' this.filename]);
%             end
            
            % advance one header
            fseek_result = fseek(this.fid,this.headBytes,'cof');
            this.fileloc = ftell(this.fid);
%             disp(['1fileloc:  ',num2str(this.fileloc) '   jmp: ',num2str(jumpBackBytes)])
            % bulk read in data
            pstring = [num2str(this.recLength) '*' this.precision];
            temp = fread( ...
                this.fid, ...
                nRecRead*this.recLength, ...
                pstring, ...
                this.footBytes+this.headBytes, ...
                this.fread_machinefmt);
            
            % jah debug
%             temp(1:10)

            % back up one header (fread skip option skips footer and next header)
            % in order to get back to record boundary, this must be done
            this.fileloc = ftell(this.fid);
            jumpBackBytes = this.fileloc - this.pageStartByte;
            fseek(this.fid,-this.headBytes-jumpBackBytes,'cof');
            this.fileloc = ftell(this.fid);
%             disp(['fileloc:  ',num2str(this.fileloc) '   jmp: ',num2str(jumpBackBytes)])
            
            if(numel(temp) == this.recLength*nRecRead)
                this.pageMemory(1:nRecRead,:) = ...
                    reshape(temp, ...
                        this.recLength, ...
                        nRecRead).';
                this.pageMemory((nRecRead+1):end,:) = 0;
            else
                nValid = numel(temp);
                nRemainingRecordsRead = nValid / this.recLength;
                this.pageMemory(1:nRemainingRecordsRead) = ...
                    reshape(temp, ...
                        this.recLength, ...
                        nRemainingRecordsRead);
                this.pageMemory((nRemainingRecordsRead+1):end,:) = 0;
            end
            
            if(direction == 0)
                this.remRecords = this.remRecords-this.pageSize*floor((1-this.pageOverlap)*this.recLength)/this.recLength;
            else
                this.remRecords = this.remRecords+this.pageSize*floor((1-this.pageOverlap)*this.recLength)/this.recLength;
            end
        end
        
        function data = getPage(this,fcn,varargin)
            if(exist('fcn'))
                if(strcmp(class(fcn),'function_handle'))
                    if(isempty(varargin))
                        data = fcn(this.pageMemory);
                    else
                        data = fcn(this.pageMemory,varargin);
                    end
                elseif(strcmp(fcn,'complex'))
                    if(size(this.pageMemory,1) == 1)
                        data = this.pageMemory(1:2:end) + 1i*this.pageMemory(2:2:end);
                    else
                        data = this.pageMemory(:,1:2:end) + 1i*this.pageMemory(:,2:2:end);
                    end
                end
            else
                data = this.pageMemory;
            end
        end
        
        function scrollData(this,fcn)
            no_fcn = ~exist('fcn');
            
            resetPage(this);
            
            figure;
            nFullPages = floor(this.nRecords/this.pageSize);
            if(no_fcn)
                hi = imagesc(getPage(this));
            else
                hi = imagesc(getPage(this,fcn));
            end
            
            ht = title({'Press any key to start.'});
            pause;
            stop = 50;
%             for k = 1 : (nFullPages-1);
            for k = 1 : stop;
                this.nextPage();
                if(no_fcn)
                    set(hi,'CData',getPage(this));
                else
                    set(hi,'CData',getPage(this,fcn));
                end
                drawnow;
            end
        end
        
        function singlePlayback(this,fcn,ylims,refidx)
            no_fcn = ~exist('fcn');
            resetPage(this);
            
            nFullPages = floor(this.nRecords/this.pageSize);
            nFullPages = 50;
            for kpage = 1 : (nFullPages-1);
                if(no_fcn)
                    dat = getPage(this);
                else
                    dat = getPage(this,fcn);
                end
                    
                if(kpage == 1)
                    figure;
                    hp = plot(dat(1,:));
                    ht = title({'Press any key to start.'},'Color',0.5*[1 1 1]);
                    if(exist('ylims'))
                        ylim(ylims);
                    end
                    xlabel('Samples');
                    ylabel('Magnitude');
                    set(gcf,'Color','k');
                    set(gca,'Color','k');
                    set(gca,'XColor',0.5*[1 1 1],'YColor',0.5*[1 1 1]);
                    set(hp,'Color','y','LineWidth',2);
                    grid on;
                    colorbar;
                    drawnow;
                    pause;
                end
                try
                    for krow = 1 : this.pageSize
                        set(hp,'YData',dat(krow,:));
                        if(~exist('refidx'))
                            set(ht,'String',['Frame ' num2str(this.pageStartRec-1+krow) ]);
                        else
                            set(ht,'String',['Reference ' num2str(dat(krow,refidx),'%3.4f') ]);
                        end
                        drawnow;
                    end
                catch
                    this.nextPage();
                end
            end
        end
        
        function reversePage(this)
            this.pageMemory = this.pageMemory(end:-1:1,:);
        end
        
        function plotimg(this,fcn,update_hi)
            if(exist('fcn'))
                if(~isa(fcn,'function_handle'))
                    error('argument must be a function handle! "@fcn" not "fcn"');
                end
                dat = fcn(this.pageMemory);
                han.fcn = fcn;
            else
                dat = this.pageMemory;
            end
            if(~exist('update_hi'))
                hf=figure;
                hi=imagesc(dat);
                colorbar;
                title({this.filename,['Frame ' num2str(this.pageStartRec) ' to ' num2str(this.pageStopRec)]});
                set(hi,'ButtonDownFcn',@cb_plot_single);
                han.hf = hf;
                han.hi = hi;
                han.this = this;
                set(hf,'UserData',han);
                set(hf,'KeyPressFcn',@cb_plot_keypress);
            else
                set(update_hi,'CData',dat);
                title({this.filename,['Frame ' num2str(this.pageStartRec) ' to ' num2str(this.pageStopRec)]});
            end
        end
        function plot(this,idxy,idxx)
            figure;
            if(exist('idxx'))
                plot(this.pageMemory(:,idxx),this.pageMemory(:,idxy));
            else
                plot(this.pageMemory(:,idxy));
            end
        end
        
        function stft(this,N)
            
%             addpath('
            
            figure
            for i = 1:N
%                 x = WaveFftObj(
            end
            
        end
    end
end