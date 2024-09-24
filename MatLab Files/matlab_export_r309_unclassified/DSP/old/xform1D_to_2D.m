function A = xform1D_to_2D(vectorData,ncolumns)
% Help: A = xform1D_to_2D(vectorData,ncolumns)
%      vectorData = row vector of samples
%      ncolumns = # of "fast time" samples
    if(size(vectorData,1) > size(vectorData,2))
        vectorData = vectorData.';
    end
    L = length(vectorData); % length of input
    nRowsF = L / ncolumns;
    nRows  = floor(nRowsF);
    Lout = nRows * ncolumns;
    
    A = reshape(vectorData(1:Lout),ncolumns,nRows).';
    
end