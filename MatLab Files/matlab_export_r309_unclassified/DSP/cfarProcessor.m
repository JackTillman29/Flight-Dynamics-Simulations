function varargout = cfarProcessor(x,alpha,nTrainingCells,nGuardCells,varargin)
% Implementation of a simple Cell Averaging (CA) Constant False Alarm Rate
% (CFAR) detector.
% Syntax:
%     Inputs:
%     varargout = cfarProcessor(x,alpha,nTrainingCells,nGuardCells,varargin)
%         ['dim', dim]:
%                dimension along which to apply the cfar threshold (in cases
%                where x is 2-D)
%         ['get_kernel', 1]:
%                this function will run until the CFAR kernel is generated and
%                will return the kernel.
%     Outputs:
%         [y] = cfarProcessor(...)
%         [y, threshold] = cfarProcessor(...)
%         [y, threshold, detections] = cfarProcessor(...)
%_____________________
% nTrainingCells and nGuardCells should be integers that sum up to the
% number of cells used in the CFAR on ONE SIDE of the cell-under-test
% DO NOT NEED TO ZERO-PAD INPUT SIGNAL.
%--------------------------------------------------------------------------
% alpha maps to linear factor relative to CFAR training mask output
%            so, for example, 13dB above average level => alpha = 10^(13/10) = 20.0
% nTrainingCells: can be scalar or a 2-element array
%                 first element is the left side of the CFAR, second
%                 element is the right side of the CFAR.
%            ex.  nTrainingCells = [2 4], row space, column space.
%          [ T T      GUARDS CUT GUARDS      T T T T ]
%            2 training cells left of CUT    4 training cells right of CUT
% nGuardCells: same format as nTrainingCells
%
% author: Jeff Hole (Booz Allen Hamilton)
% 7-1-2018

% SOME BASIC MATH FOR CFAR KERNEL FOOTPRINTS:
%       Tx = ntrain(1) (all one-sided)
%       Ty = ntrain(2)
%       Gx = nguard(1)
%       Gy = nguard(2)
%       total cfar footprint = 
%           (2*(Tx+Gx) + 1) * (2*(Ty+Gy) + 1)
%                   ** the "+1" is for the CUT
%           = (2Tx + 2Gx + 1) * (2Ty + 2Gy + 1)
%              [EXPAND]
%           = 4TxTy + 4TxGy + 2Tx +
%             4GxTy + 4GxGy + 2Gx +
%             2Ty   + 2Gy   + 1
%        total guard + CUT footprint = 
%           (2*Gx + 1) * (2*Gy + 1)
%              [EXPAND]
%           = 4GxGy + 2Gx + 2Gy + 1
%        AVERAGING FACTOR = total cfar footprint - "guard+CUT" footprint
%           = 4TxTy + 4TxGy           + 2Tx +
%             4GxTy + (4GxGy - 4GxGy) + (2Gx - 2Gx) +
%             2Ty   + (2Gy-2Gy)       + (1 - 1)
%           = 4TxTy + 4TxGy + 2Tx +
%             4GxTy + 0     + 0 +
%             2Ty   + 0     + 0

% NEED TO HANDLE EDGES EXPLICITLY
%   C = cell-under-test
%   G = guard cell
%   T = training cell
%       edge
%       -----------------------------
%       | T | G | G | G | T | T |  
%       ---------***---------------
% edge  | T | G |*C*| G | T | T |  
%       ---------***---------------
%       | T | G | G | G | T | T |     ...
%       ---------------------------
%       | T | T | T | T | T | T |  
%       ---------------------------
%       | T | T | T | T | T | T |  
%       ---------------------------
%       |   |   |   |   |   |   |  
%       |             .
%       |             .
%       |             .
%    2 KINDS OF BOUNDARIES:
%        1) HYPOTHETICAL CELLS OUTSIDE THE BORDERS ARE SET TO ZERO.
%        2) CONNECT BORDERS TO CREATE PERIODIC BOUNDARIES (TOP-BOTTOM, 
%           LEFT-RIGHT forms a torus)
%        3) NO CELLS OUTSIDE THE BORDER (DEFAULT)
%			Adaptively calculate the number of active averaging cells at the borders
%			by computing the conv() between the 
%    FOR CFAR, WE SHOULD BE ABLE TO CAPTURE ANY KIND OF BOUNDARIES WE WANT.
%    ** FOR SCANNING RADARS, A PERIODIC BOUNDARY MIGHT BE BEST (2ND KIND)
%    ** FOR 2D RDM, BEST TO JUST USE THE ACTUAL CELLS (3RD KIND)

% TO IMPLEMENT 3RD KIND:
%   ON EDGES, NEED TO ONLY SUM OVER CELLS WITHIN THE BOUNDARY AND AVERAGE
%   BY THE NUMBER OF ACTUAL CELLS (DO NOT ASSUME CELLS BEYOND THE
%   BOUNDARIES ARE EQUAL TO ZERO)
%     -- COULD RUN MATLAB conv2, BUT THEN SCALE THE AVERAGING FACTOR FOR
%        EACH CUT USING AN "AVERAGING FACTOR MASK"
%        
%        EXAMPLE
%            ntrain = [1 1] % one-sided
%            nguard = [0 0]
%        total cfar footprint = (2*(1+0)+1) * (2*(1+0)+1) = 16
%                             = (2+1) * (2+1) = 9
%        ----------------------------------
%        |   |   |   |   |   |   |   |
%        ----------------------------------
%        |   |   |   |   |   |   |   |
%        ----------------------------------
%        |   |   |   |   |   |   |   |
%        ----------------------------------
%        |   |   |   |   |   |   |   |
%        ----------------------------------
%        |   |   |   |   |   |   |   |
%        ----------------------------------
%        |   |   |   |   |   |   |   |



% N = length(x);
% nCells = 2*nTrainingCells;
% alpha = nCells .* (Pfa^(-1/nCells) - 1);
% y          = zeros(N,1);
% threshold  = zeros(N,1);
% detections = zeros(N,1);


get_kernel = 0;
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case {'dim','dimension','d'}
            dim = varargin{k+1};
        case 'get_kernel'
            get_kernel = varargin{k+1};
    end
end
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
    
    if(get_kernel)
        varargout{1} = cfar_kernel;
        return
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
        % TODO: ensure the these implement the default (no cells outside signal, 
		%       i.e. do not assume zeros outside, but adaptively calculate # of 
		%       averaging cells on the border
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
else % Checked Hole/Sawmiller 9/23/2019

    iCUT = nTrainingCells+nGuardCells+1;
    iGuard1 = iCUT-nGuardCells; % [row space,column space]
    iGuard2 = iCUT+nGuardCells; % [row space,column space]
    cfar_kernel = ones(2*nTrainingCells+2*nGuardCells+1);
    cfar_kernel(iGuard1(1):iGuard2(1),iGuard1(2):iGuard2(2)) = 0;
%     cfar_kernel = cfar_kernel ./ sum(cfar_kernel(:)); % normalize

    if(get_kernel)
        varargout{1} = cfar_kernel;
        return
    end
    
    % this line effectively counts the number of training cells active for
    % each cell-under-test. Dividing the output of the convolution of the
    % input with the cfar footprint effectively computes the average noise
    % estimate, and handles the edges where otherwise the threshold would
    % drop yielding false detections on the edges.
    num_training_cells = conv2(ones(size(x)),cfar_kernel,'same');
    
    cfar = conv2(x,flipup(fliplr(cfar_kernel)),'same') ./ num_training_cells;
    
    if(debug)
        figure; imagesc(cfar_kernel); colorbar; title('CFAR Kernel')
        axis equal; axis tight;
%         imagesc(num_training_cells);
    end

end




threshold  = alpha * cfar;
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