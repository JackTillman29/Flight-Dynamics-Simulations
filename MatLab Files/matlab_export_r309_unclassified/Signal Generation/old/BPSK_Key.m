function [key,length,keyu,bitkey] = BPSK_Key(nChips,fbr)
    % nChips = 2^n-1
    key = zeros(nChips+1,1);
    nStages = log10(nChips+1) / log10(2);
    disp([num2str(nStages) ' stages needed for this sequence length']);
    stages = uint8(ones(1,round(nStages)));
    max_found = 0;
    
    if(max(fbr) > (nStages ))
        error(['Only ' num2str(nStages) ' needed for this sequence length!']);
    end
    key(1) = uint64(uintArray2Int(stages));
    %disp([key(1) uint64(stages)]);
    for k = 2 : (nChips + 2)
        %feedback = mod(stages(fbr(1)) + stages(fbr(2)),2);
        feedback = mod(sum(stages(fbr)),2);
        
        stages = [ feedback stages(1:(end-1))];
        val = uint64(uintArray2Int(stages));
        if( (val == key(2)) && (max_found == 0) )
            length = k - 2;
            max_found = 1;
            max_found_j = k;
        end
        key(k) = val;
        %disp([key(k) uint64(stages)]);
    end
    keyu = key;
    try
        key = key(2:(max_found_j-1));
    catch
        key_reversed = key(end:-1:1);
        [a,ai] = find(key_reversed == key_reversed(1));
        length = a(2) - 1;
        key = key_reversed([(a(2)-1):-1:a(1)]);
    end
    bitkey = bitand(key,1);
end