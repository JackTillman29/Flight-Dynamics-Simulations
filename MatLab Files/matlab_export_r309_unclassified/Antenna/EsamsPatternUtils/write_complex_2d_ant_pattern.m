function esamsTable = write_complex_2d_ant_pattern(Az,El,data,filename)

%use this to write complex antenna patterns in esams format
%for example, a complex element pattern which may be read 
%in by the array code (which expects angles in radians)
%
% this function is similar to write_2d_ant_pattern.m (which is magnitude
% only)
%
% Az - vector of az points for data array (- to +)
% El - vector of el points for data array (- to +)
% data - complex antenna pattern, with dimension length(El) x length(Az)

data(isnan(data)) = 0.01;

%force these to be row vectors
Az = Az(:).';
El = El(:).';

tempData1 = [Az;                real(data);];
tempData2 = [zeros(size(Az));   imag(data);];

esamsTable(:,1) = [length(Az) length(El) El tempData1(:).']';
esamsTable(:,2) = [0 0 zeros(size(El)) tempData2(:).']';

fid = fopen(filename,'w+');
fprintf(fid,'%f\t%f\n',esamsTable.');
fclose(fid);

end