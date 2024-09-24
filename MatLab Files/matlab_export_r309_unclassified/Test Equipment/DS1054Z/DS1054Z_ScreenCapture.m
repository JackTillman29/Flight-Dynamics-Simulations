function data_object = DS1054Z_ScreenCapture(obj,flag_invert,outFilename)
%% This function requires that the SCPI I/F be "closed"
fclose(obj);
warning('off','instrument:fread:unsuccessfulRead');

%% Get current state
fopen(obj)

fprintf(obj,':STOR:IMAG:TYPE PNG');
if(flag_invert)
    fprintf(obj,':STOR:IMAG:INVERT ON');
else
    fprintf(obj,':STOR:IMAG:INVERT OFF');
end

% ======================= channel configuration
fprintf(obj,':DISP:DATA?')
dat = zeros(1,1152066,'uint8');
idx1 = 1;
for k = 1 : 4
    
    rdat = char(fread(obj));
    idx2 = idx1 - 1 + 250012;
    dat(idx1:idx2) = rdat;
    idx1 = idx2 + 1;
    if(k == 1)
        disp(['Rx Header: ' char(dat(1:11))]);
        disp(['Image Type: ' char(dat(12:15))]);
    end
end
idx2 = idx1 - 1 + 152018;
rdat = char(fread(obj));
dat(idx1:idx2) = rdat;

fido = fopen(outFilename,'wb');
fwrite(fido,dat(12:end));
fclose(fido);



fclose(obj);
end