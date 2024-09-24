function [x y values] = esamsToMatlabTable(esamsArray)
numX = esamsArray(1);
numY = esamsArray(2);

x = esamsArray(numY+3:numY+1:numX+numY+2+numX*numY);
y = esamsArray(3:numY+2);

k = 1;
for j = numY + 3:numY+1:numX+numY+2+numX*numY
    for l = 1:numY
        values(k,l) =esamsArray(j+l);
    end
    k = k+1;
end