function C = convmatrix(kernel,rxSignal,mode)

if(~exist('mode','var'))
    warning('Returning full convolution matrix. You may pass in ''full'' or ''trimmed''');
    mode = 'full';
end

M = length(rxSignal);
N = length(kernel);
L = M + N - 1;

% Create convolution matrix (Toeplitz array)
C = zeros(L,M,class(kernel));
for k = 1 : M
    C(k:(k+N-1),k) = kernel(:)';
end

switch(lower(mode))
    case 'full'
        % do nothing
    case 'trimmed'
        % return only the section of the array that has a full set of
        % coefficients.
        row_start = floor(length(kernel)/2);
        row_stop  = L - row_start;
        C = C(row_start:row_stop,:);
    otherwise
        error(['Unknown mode passed: ' mode]);
end

end