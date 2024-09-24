function outputMat = embedMatIntoMat(insertMat,targetMat,targetIdx)

% centerIdx = round(size(insertMat)/2);

N = size(targetMat);

wid1 = size(insertMat,1);
wid2 = size(insertMat,2);

newidx1 = targetIdx(1) + [0:(wid1-1)] - floor(wid1/2);
newidx2 = targetIdx(2) + [0:(wid2-1)] - floor(wid2/2);

newidx1 = newidx1((newidx1 > 0) & (newidx1 <= N(1)));
newidx2 = newidx2((newidx2 > 0) & (newidx2 <= N(2)));

outputMat = targetMat;

if ~isempty(newidx1) & ~isempty(newidx2) 
    small_idx1 = newidx1 - targetIdx(1) + 1 + floor(wid1/2);
    small_idx2 = newidx2 - targetIdx(2) + 1 + floor(wid1/2);
    outputMat(newidx1,newidx2) = insertMat(small_idx1,small_idx2);
end


end