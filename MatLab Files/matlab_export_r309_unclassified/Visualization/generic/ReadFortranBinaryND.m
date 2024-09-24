function data = ReadFortranBinaryND(inputfile,dims,varargin)
% inputFile - binary file to read in 
% dims - RxCxD 
% 'type', can be 'real' (default) or 'complex'
% 'size', can be 'single' (default) or 'double'
% ex: data = ReadFortranBinaryND('infile.bin',[2 3 4],'size','double','type','complex')
% =============== Requires the write format shown below ================
% OPEN(231, FILE = FileNameOut, ACCESS = "STREAM", FORM = "UNFORMATTED")
% WRITE(231) IQData%ClutterSignals
    
    nPairs = length(varargin);
    if(mod(nPairs,2) ~= 0)
        error('incorrect input pairs');
    end
    
    readType = 'single'; % default
    nBytesPerValue = 4; % default
    isComplex = 0; % default
    
    for k = 1:2:nPairs
        switch lower(varargin{k})
            case 'type'
                switch lower(varargin{k+1})
                    case 'real'
                        isComplex = 0;
                    case 'complex'
                        isComplex = 1;
                    otherwise
                        error(['Unknown parameter given: ' varargin{k} ' -> ' varargin{k+1}]);
                end
            case 'size'
                switch lower(varargin{k+1})
                    case 'single'
                        nBytesPerValue = 4;
                        readType='single';
                    case 'double'
                        nBytesPerValue = 8;
                        readType='double';
                    otherwise
                        error(['Unknown parameter given: ' varargin{k} ' -> ' varargin{k+1}]);
                end
            otherwise
                error(['Unknown parameter given: ' varargin{1}]);
        end
    end
        
    fid = fopen(inputfile,'rb');
    if(fid == -1)
        error(['could not locate file: ' inputfile]);
    end
    fseek(fid,0,'eof');
    nBytesTotal = ftell(fid);
    nValuesTotal = nBytesTotal / nBytesPerValue;
    
    if(mod(nValuesTotal,1) ~= 0)
        error(['Could not interpret # of writes!']);
    end
    
    frewind(fid);
    
    % reads in all data if real, real samples if complex
    data = fread(fid,nValuesTotal,readType);
    
    if(isComplex)
        data = data(1:2:end) + 1i * data(2:2:end);
    end
    
    
    fclose(fid);
    
    if(~exist('dims'))
        dims = [nValuesTotal 1 1];
    end
    data = reshape(data,dims(1),dims(2),dims(3));
    if(strcmp(readType,'single'))
        data = single(data);
    end
end