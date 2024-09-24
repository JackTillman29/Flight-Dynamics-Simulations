function R = smpcormtx(x,N)
    % Function expects temporal data to be along 2nd dimension.
    % Sample sets are along 1st dimension.
    % Data format: (elements x time samples)
    if(~exist('N'))
        N = size(x,2);
    end
    
    R = (1/N) * (x(:,1:N) * x(:,1:N)');
    
end