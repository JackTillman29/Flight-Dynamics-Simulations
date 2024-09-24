function y = rms(x)
    % perform rms calculation along the longest dimension
    dim = 1;
    if(size(x,2) > size(x,1))
        dim = 2;
    end
    
    y = sqrt(mean(abs(x).^2,dim));
end