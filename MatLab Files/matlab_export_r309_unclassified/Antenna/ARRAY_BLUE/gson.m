function E = gson(orig_basis,null_space)
    % nullspace must be column vector
    nCols = size(orig_basis,2);
    
    % normalize input basis
    for k = 1 : nCols
        orig_basis(:,k) = orig_basis(:,k) ./ norm(orig_basis(:,k));
    end
    
    
    if(exist('null_space','var'))
        % project desired nullspace onto basis
        proj = orig_basis'*null_space;
        [vmax,imax]=max(abs(proj));

        % replace basis
        rmBasis = orig_basis(:,imax);
        aPrime = null_space./norm(null_space);
        orig_basis(:,imax) = orig_basis(:,1); % move first basis here so that
        % null space vector passed in can go in column 1
        
        orig_basis(:,1) = aPrime;
        
    end
    
    % now perform gram-schmidt orhtogonalization
    E = zeros(size(orig_basis));

    for k = 1 : nCols 
        %E(:,k)
        temp = 0;
        for ki = 1 : (k-1)
            temp = temp + E(:,ki)' * orig_basis(:,k) * E(:,ki);
        end
        E(:,k) = orig_basis(:,k) - temp;
        E(:,k) = E(:,k) ./ norm(E(:,k));
    end 
    
    if(exist('null_space','var'))
        E = E(:,2:end);
    end
    
end