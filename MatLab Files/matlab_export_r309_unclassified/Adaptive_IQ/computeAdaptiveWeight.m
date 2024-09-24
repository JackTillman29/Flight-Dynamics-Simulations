function w = computeAdaptiveWeight(s,U,V,mode,modeval)
    % s is quiescent steering vector
    % U is vector of eigenvalues of correlation matrix, R
    % V is eigenvector set associated with U
    % Umin is minimum eigenvalue to use
    
    switch(lower(mode))
        case 'direct'
            N = modeval;
            Umin = 0; % ensure Umin is less than U(N)
        case 'mineig'
            % How many terms will be have?
            N = sum(diag(U) > Umin);
        otherwise
            error('Unknown mode passed.');
    end
    
   
    
    
    
    w = s;
    for k = 1 : N
        w = w - ((U(k)-Umin)./U(k)) * ((V(:,k)' * s') * V(:,k))';
    end
    
    % normalize the output
    %w = w ./ abs(w);
end

%s - (V(:,1)' * s') * V(:,1)'