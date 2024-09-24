function y = movavg(x,n)
len = size(x,1);
nVectors = size(x,2);
numAves = n;
y = zeros(size(x,1),size(x,2));
y(1:numAves-1) = x(1:numAves-1);
for ki = 1:nVectors
    for k = numAves:len %outer loop
        y(k,ki) = (x(k,ki) + sum(x(k-[1:numAves-1],ki)))./numAves;
    end
end