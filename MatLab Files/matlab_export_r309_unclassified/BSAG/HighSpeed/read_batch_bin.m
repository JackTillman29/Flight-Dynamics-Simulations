% fwrite(fidbw,[size(imData,2) sum(imData,1)./size(imData,1) sum(imData_TG,1)./size(imData_TG,1)],'single');
fid = fopen('c:\temp\avgdata1.bin','rb');
nWrites = 0;
nVals = [];
fseek_stat = 0;
while(~fseek_stat)
    n = fread(fid,1,'single');
    nVals = [nVals n];
    if(~isempty(n))
        fseek_stat = fseek(fid,4*(2*n),'cof');
    else
        fseek_stat = 1;
    end
    if(fseek_stat==0)
        nWrites = nWrites + 1;
    end
    disp([size(nVals) nWrites]);
end

frewind(fid);



fclose(fid);