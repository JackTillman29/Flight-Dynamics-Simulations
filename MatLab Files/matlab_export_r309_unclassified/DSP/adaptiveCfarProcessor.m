function varargout = adaptiveCfarProcessor(x,numStdDev,nTrainingCells,nGuardCells,varargin)
% varargout = adaptiveCfarProcessor(x,numStdDev,nTrainingCells,nGuardCells,[dim])
%     [dim]: dimension along which to apply the cfar threshold (in cases
%            where x is 2-D)
%  ** OUTPUT threshold is the combined CFAR or "full" CFAR which uses
%  ** both sides of the cell-under-test (CUT) together and and is
%  ** based on the equation MEAN + numStdDev*STDDEV
% OUTPUT OPTIONS:
%     [y] = cfarProcessor(...)
%     [y, threshold] = cfarProcessor(...)
%     [y, threshold, detections] = cfarProcessor(...)
%     [y, threshold, detections, side_thresholds_LR] = adaptiveCfarProcessor(...)
%         side_thresholds_LR is a cell array containing the single-sided
%         thresholds (left/early and right/late sides)

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
%
% author: Jeff Hole (Booz Allen Hamilton)
% 7-1-2018

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
        cfar_kernel_R = [ ...
            zeros(1,nTrainingCells(1)) ...
            zeros(1,nGuardCells(1)) ...
            0 ... % CUT
            zeros(1,nGuardCells(2)) ...
            ones(1,nTrainingCells(2))]/nTrainingCells(2);
        cfar_kernel_L = [ ...
            ones(1,nTrainingCells(1)) ...
            zeros(1,nGuardCells(1)) ...
            0 ... % CUT
            zeros(1,nGuardCells(2)) ...
            zeros(1,nTrainingCells(2))]/nTrainingCells(1);
        cfar_kernel_LR = cfar_kernel_L + cfar_kernel_R;
    % unbalanced
    else
        nTrainMax = max(nTrainingCells);
        nGuardMax = max(nGuardCells);
        nMaxCells = nTrainMax + nGuardMax;
        % pad ends with zeros so that CUT is always in the middle (+/- 1
        % sample) of the kernel
        cfar_kernel_R = [ ...
            zeros(1,nMaxCells-nTrainingCells(1)-nGuardCells(1)) ...
            zeros(1,nTrainingCells(1)) ...
            zeros(1,nGuardCells(1)) ...
            0 ... % CUT
            zeros(1,nGuardCells(2)) ...
            ones(1,nTrainingCells(2)) ...
            zeros(1,nMaxCells-nTrainingCells(2)-nGuardCells(2)) ]/(sum(nTrainingCells));
        cfar_kernel_L = [ ...
            zeros(1,nMaxCells-nTrainingCells(1)-nGuardCells(1)) ...
            ones(1,nTrainingCells(1)) ...
            zeros(1,nGuardCells(1)) ...
            0 ... % CUT
            zeros(1,nGuardCells(2)) ...
            zeros(1,nTrainingCells(2)) ...
            zeros(1,nMaxCells-nTrainingCells(2)-nGuardCells(2)) ]/(sum(nTrainingCells));
        cfar_kernel_LR = cfar_kernel_L + cfar_kernel_R;
    end
    
    ncfar = length(cfar_kernel_R);
    if( any(size(x)==1) )
        
        % perform conv with a moving avg filter effectively computes the
        % MOVING AVG within the active taps of the moving avg filter
        %    zero-padded movavg filters provides the proper offset rel to
        %    CUT for measuring the mean and std of the training/averaging
        %    cells
        mean_L  = conv(x,cfar_kernel_L(end:-1:1),'same');
        mean_R  = conv(x,cfar_kernel_R(end:-1:1),'same');
        mean_LR = conv(x,cfar_kernel_LR(end:-1:1),'same');
        
        % E[x] = mean(x), E[x^2] = mean(x^2), E[x^2] - E[x]^2 = var(x)
        Ex2_L  = conv(x.^2,cfar_kernel_L(end:-1:1),'same');
        Ex2_R  = conv(x.^2,cfar_kernel_R(end:-1:1),'same');
        Ex2_LR = conv(x.^2,cfar_kernel_LR(end:-1:1),'same');
        
        var_L  = Ex2_L - mean_L.^2;
        var_R  = Ex2_R - mean_R.^2;
        var_LR = Ex2_LR - mean_LR.^2;
        
        std_L  = sqrt(var_L);
        std_R  = sqrt(var_R);
        std_LR = sqrt(var_R);
        
        
        
        
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
        
        % perform the conv of the cfar kernel on a signal matrix across a
        % given dimension
        % "dim" is a varargin
        if(dim == 1)
            
        % cfar across columns
        elseif(dim == 2)
            
            xpad = [zeros(size(x,1),floor(ncfar/2)) x zeros(size(x,1),floor(ncfar/2))]
            
            mean_L = filter(cfar_kernel_L(end:-1:1),1,xpad,[],dim); mean_L = mean_L(:,ncfar:end);
            mean_R = filter(cfar_kernel_R(end:-1:1),1,xpad,[],dim); mean_R = mean_R(:,ncfar:end);
            
            Ex2_L = filter(cfar_kernel_L(end:-1:1),1,xpad.^2,[],dim); Ex2_L = Ex2_L(:,ncfar:end);
            Ex2_R = filter(cfar_kernel_R(end:-1:1),1,xpad.^2,[],dim); Ex2_R = Ex2_R(:,ncfar:end);
            
            var_L = Ex2_L - mean_L.^2;
            var_R = Ex2_R - mean_R.^2;

            std_L = sqrt(var_L);
            std_R = sqrt(var_R);
            
            
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




threshold_L  = mean_L  + numStdDev * std_L;
threshold_R  = mean_R  + numStdDev * std_R;
threshold_LR = mean_LR + numStdDev * std_LR;
% mean_LR = min(mean_L,mean_R);
% threshold = 

% choose the lower threshold
threshold  = min(threshold_L,threshold_R);
detections = (x>threshold).*1; % 0 or 1
y          = detections.*x;     % pass original signal value through

% detections_L = (x>threshold_L).*1;  % 0 or 1
% detections_R = (x>threshold_R).*1;  % 0 or 1




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
elseif(nargout == 4)
    varargout{1} = y;
    varargout{2} = threshold;
    varargout{3} = detections;
    varargout{4} = {threshold_L, threshold_R, threshold_LR};
end


end