function data = ReadFortranBinaryDouble(inputfile)
    % Assumes all data is single precision floating point, 4 bytes per
    % float
    fid = fopen(inputfile,'rb');
    
    
    nBytesPerWrite = fread(fid,1,'uint32');
    nValuesPerWrite = nBytesPerWrite / 8.0;
    fseek(fid,0,'eof');
    totalLines = ftell(fid)/(8+nBytesPerWrite);
    frewind(fid);
       
    data = zeros(totalLines,nValuesPerWrite);

    for k = 1:totalLines
        fread(fid,1,'uint32');
        data(k,:) = fread(fid,nValuesPerWrite,'double');
        fread(fid,1,'uint32');
    end
    fclose(fid);
end