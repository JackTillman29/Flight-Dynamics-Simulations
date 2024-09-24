function varargout = cfarProcessor_log(x,alpha,nTrainingCells,nGuardCells,varargin)
% varargout = cfarProcessor(x,alpha,nTrainingCells,nGuardCells,[dim])
%     [dim]: dimension along which to apply the cfar threshold (in cases
%            where x is 2-D)
%_____________________
% nTrainingCells and nGuardCells should be integers that sum up to the
% number of cells used in the CFAR on ONE SIDE of the cell-under-test
% DO NOT NEED TO ZERO-PAD INPUT SIGNAL.
%--------------------------------------------------------------------------
% nTrainingCells: can be scalar or a 2-element array
%                 first element is the left side of the CFAR, second
%                 element is the right side of the CFAR.
%            ex.  nTrainingCells = [2 4]
%          [ T T      GUARDS CUT GUARDS      T T T T ]
%            2 training cells left of CUT    4 training cells right of CUT
% nGuardCells: same format as nTrainingCells


% N = length(x);
% nCells = 2*nTrainingCells;
% alpha = nCells .* (Pfa^(-1/nCells) - 1);
% y          = zeros(N,1);
% threshold  = zeros(N,1);
% detections = zeros(N,1);

if(length(varargin) > 0)
    dim = varargin{1};
elseif(any(size(x)==1))
    dim = 1;
else
    dim = [];
end

% if single element
if(all(size(nTrainingCells) == 1))
    % duplicate to make two elements
    nTrainingCells = [nTrainingCells nTrainingCells];
end

% if single element
if(all(size(nGuardCells) == 1))
    % duplicate to make two elements
    nGuardCells = [nGuardCells nGuardCells];
end

k1CUT = nTrainingCells(1) + nGuardCells(1) + 1;
k2CUT = nTrainingCells(2) + nGuardCells(2) + 1;

debug = 0;

% 1-D CFAR
if( ~isempty(dim) )
    
    % balanced
    if(k1CUT == k2CUT)
        cfar_kernel = [ ...
            ones(1,nTrainingCells(1)) ...
            zeros(1,nGuardCells(1)) ...
            0 ... % CUT
            zeros(1,nGuardCells(2)) ...
            ones(1,nTrainingCells(2))]/(sum(nTrainingCells));
    % unbalanced
    else
        nTrainMax = max(nTrainingCells);
        nGuardMax = max(nGuardCells);
        nMaxCells = nTrainMax + nGuardMax;
        % pad ends with zeros so that CUT is always in the middle (+/- 1
        % sample) of the kernel
        cfar_kernel = [ ...
            zeros(1,nMaxCells-nTrainingCells(1)-nGuardCells(1)) ...
            ones(1,nTrainingCells(1)) ...
            zeros(1,nGuardCells(1)) ...
            0 ... % CUT
            zeros(1,nGuardCells(2)) ...
            ones(1,nTrainingCells(2)) ...
            zeros(1,nMaxCells-nTrainingCells(2)-nGuardCells(2)) ]/(sum(nTrainingCells));
    end
    
    ncfar = length(cfar_kernel);
    if( any(size(x)==1) )
        cfar = conv(x,cfar_kernel(end:-1:1),'same');
%         % uneven CFAR
%         if(k1CUT ~= k2CUT)
%             cfar = reshape(cfar,1,[]);
%             if(k1CUT > k2CUT)
%                 kdelta = floor((k1CUT-k2CUT)/2);
%                 cfar = [zeros(1,kdelta) cfar(1:(end-kdelta))];
%             else
%                 kdelta = floor((k2CUT-k1CUT)/2);
%                 cfar = cfar(end:-1:1);
%                 cfar = [zeros(1,kdelta) cfar(1:(end-kdelta))];
%                 cfar = cfar(end:-1:1);
%                 cfar = reshape(cfar,size(x));
%             end
%         end
    else
        % work-in-progress
        
        % 
        if(dim == 1)
            
        % cfar across columns
        elseif(dim == 2)
            cfar = filter(cfar_kernel(end:-1:1),1,...
                [zeros(size(x,1),floor(ncfar/2)) x zeros(size(x,1),floor(ncfar/2))],[],dim);
            cfar = cfar(:,ncfar:end);
        end
        
        
%         cfar = filter(cfar_kernel,1,x,[],dim);
%         N = size(x,dim);
%         M = size(cfar_kernel,dim);
%         L = N + M - 1
%         
%         cfar = ifft( fft(cfar_kernel,L,dim) .* fft(x,L,dim) );
%         if(dim == 2)
%             cfar = [ ...
%                  cfar(:,(L-M+1):end) ...
%                  cfar(:,1:(L-M)) ];
%             K2 = floor(M/2);
%             if(length(K2:(L-K2)) ~= size(x,dim))
%                 cfar = cfar(:,K2:(end-K2-1)) ./ M;
%             else
%                 cfar = cfar(:,K2:(end-K2)) ./ M;
%             end
%         elseif(dim == 1)
%             cfar = [ ...
%                  cfar((L-M+1):end,:); ...
%                  cfar(1:(L-M),:) ];
%             K2 = floor(M/2);
%             if(length(K2:(L-K2)) ~= size(x,dim))
%                 cfar = cfar(K2:(end-K2-1),:) ./ M;
%             else
%                 cfar = cfar(K2:(end-K2),:) ./ M;
%             end
%         end
    end
% 2-D CFAR
else
    nCells = 2*nTrainingCells+2*nGuardCells+1;
    iCUT = nTrainingCells+nGuardCells+1;
    iGuard1 = iCUT-nGuardCells;
    iGuard2 = iCUT+nGuardCells;
    cfar_kernel = ones(nCells)./((2*nTrainingCells).^2-(2*nGuardCells).^2 - 1);
    cfar_kernel(iGuard1:iGuard2,iGuard1:iGuard2) = 0;
    cfar = conv2(x,cfar_kernel,'same');
    
    if(debug)
        figure; imagesc(cfar_kernel); colorbar; title('CFAR Kernel')
    end

end




threshold  = alpha + cfar; % ADD BECAUSE IT IS LOG INPUT SIGNAL
detections = (x>threshold).*1;  % 0 or 1
y          = detections.*x;     % pass original signal value through

% figure; imagesc(detections); colorbar;


% % EXPLICIT METHOD
% for k = 1:N
%     if(k-nTrainingCells-nGuardCells <= 0)
%         % BEGINNING OF SIGNAL ARRAY
%         cfar = ( sum(x(1:(k+nTrainingCells+nGuardCells))) - ...
%             x(k) - ...
%             sum(x(max(k-nGuardCells,1):max(k-1,1))) - ...
%             sum(x((k+1):(k+nGuardCells))) ) / (nCells-(nTrainingCells-k));
%     elseif(k+nTrainingCells+nGuardCells >= N)
%         % END OF SIGNAL ARRAY
%         cfar = ( sum(x((k-nTrainingCells-nGuardCells):end)) - x(k) ) / (nCells-(nTrainingCells-(N-k)));
%     else
%         % REST OF SIGNAL ARRAY
%         cfar = ( sum(x([-(nTrainingCells+nGuardCells):(nTrainingCells+nGuardCells)]+k)) - ... % add all cells nTrain+nGuard before and after CUT and CUT
%             x(k) - ... % subtract out CUT power
%             sum(x((k-nGuardCells):(k-1))) - ... % subtract leading guard cells
%             sum(x((k+1):(k+nGuardCells))) ) / nCells; % subtract lagging guard cells
%     end
%     
%     threshold(k) = cfar * alpha;
%     if(x(k) > threshold(k))
%         y(k) = x(k);
%         onoff(k) = 1;
%     else
%         y(k) = 0;
%         onoff(k) = 0;
%     end
% end



if(nargout == 1)
    varargout{1} = y;
elseif(nargout == 2)
    varargout{1} = y;
    varargout{2} = threshold;
elseif(nargout == 3)
    varargout{1} = y;
    varargout{2} = threshold;
    varargout{3} = detections;
end


end