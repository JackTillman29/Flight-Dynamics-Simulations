function [detX,detY] = cfarRDM_Merge(detMatrix)
    % expects output from cfarProcessor "det"
    
    % build vector of vertical (e.g. range) detections
    rngSumDets = sum(detMatrix,1);
    
    % build a list of columns which have more than 1 detection
    idxColumns = find(rngSumDets > 1);
    detX = [];
    detY = [];
    
    % loop over columns with more than 1 detection
    for k = 1 : length(idxColumns)
        fprintf('Checking column %d\n',idxColumns(k));
        iStart = 0;
        for p = 1 : size(detMatrix,1)
            
            if(detMatrix(p,idxColumns(k)) == 1)
                if(iStart == 0)
                    iStart = 1;
                    mergeRowStart = p;
                end
            else
                if(iStart == 1)
                    % end of the line
                    iStart = 0;
                    
                    detLength = p - mergeRowStart;
                    mergeRow = mergeRowStart + detLength/2-0.5;
                    detX = [detX idxColumns(k)];
                    detY = [detY mergeRow];
                end
            end
        end
        if(iStart == 1) % got to end of row and did not terminate a merge
            detLength = p - mergeRowStart;
            mergeRow = mergeRowStart + detLength/2;
            detX = [detX idxColumns(k)];
            detY = [detY mergeRow];
        end
    end
    
    
    

end

