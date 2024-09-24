function data = ReadFortranBinary(inputfile,varargin)
    % Assumes all data is single precision floating point, 4 bytes per
    % float
    fid = fopen(inputfile,'rb');
    
    nBytesPerWrite = fread(fid,1,'uint32')
    nValuesPerWrite = nBytesPerWrite / 4.0
    fseek(fid,0,'eof');
    totalLines = ftell(fid)/(8+nBytesPerWrite);
    frewind(fid);
    if(nargin == 1 || strcmp(varargin{1},'single'))
        data = single(zeros(totalLines,nValuesPerWrite));
    elseif(strcmp(varargin{1},'double'))
        data = zeros(totalLines,nValuesPerWrite);
    else
        error(['Argument 2 of ReadFortranBinary is invalid! Must be either single, double, or omitted.']);
    end
    for k = 1:totalLines
        fread(fid,1,'uint32');
        data(k,:) = fread(fid,nValuesPerWrite,'single');
        fread(fid,1,'uint32');
    end
    fclose(fid);
end