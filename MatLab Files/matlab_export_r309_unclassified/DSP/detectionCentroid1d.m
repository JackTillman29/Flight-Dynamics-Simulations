function varargout = detectionCentroid1d(dets,detpwr)
% Computes the centroids (weighted and/or unweighted by detection power in
% detected cells) of connected regions
% Inputs:
%   dets:              1xN or Nx1 binary data (common output from CFAR)
%   detpwr [optional]: 1xN or Nx1 power data used to weight the centroid of
%                      cells in binary data "dets"
%
% data:   0 0 0 0 1 1 1 0 0 1 1 0 1 1 1 1 1
% region:           1        2        3
%
% author: Jeff Hole (Booz Allen Hamilton), 2019-12

if(~exist('detpwr'))
    % make "detpwr" equal to all ones for equal-weighted averaging of centroid
    detpwr = 0*dets + 1;
end

ilabel = 0;
icombined = 0; % not a component of the core algorithm
% Hoshen-Kopelman Algorithm
label = 0.*dets;
N = length(dets);
for ix = 1:N(1)
    if(dets(ix) == 1)

        if(ix == 1)
            left = 0;
        else
            left = dets(ix-1);
        end

        if( (left == 0) )
            % if not left and not above
            %   THEN NEW CLUSTER FOUND
            ilabel = ilabel + 1;
            label(ix) = ilabel;
        else
            icombined = icombined + 1;
            % if left and above are already labeled, then the CUT
            % connects these two regions, so need to make sure they are
            % all labeled the same (prioritizing label assigned to the
            % LEFT cell)
            label_left = label(ix-1);
            label(ix) = label_left;
        end

    else

    end
end

maxRegionNumber = max(label(:));

% compute centroids of each region
%    to do this, compute the mean of the region's pixel locations

ivalid_region_counter = 0;
for iregion = 1:maxRegionNumber
    
    % binary matrix with 1 where only this region exists
    % (all other regions are zeroed out)
    this_region = 1*(label == iregion);
    
    this_ind = find(label == iregion);
    
    % WORKS WITHOUT THIS IF(~ISEMPTY(this_ind))
    if(~isempty(this_ind))
        
        ivalid_region_counter = ivalid_region_counter + 1;
        
        this_pwr = detpwr(this_ind);
        
        % center-of-mass analogy:
        %    summation(i=1,N, (xi-CM)*mi) = 0
        %      solve for CM ("Center-of-Mass"):
        %           (x1*m1 - CM*m1) + ... + (xN*mN - CM*mN) = 0
        %           sum(i=1,N, xi*mi) = CM * sum(i=1,N, mi)
        %           CM = sum(i=1,N, xi*mi) / sum(i=1,N, mi)
        
        % if mi = 1 for all i, then CoM equation reduces to simply
        % computing the mean of the indices
        centroid_ind = mean(this_ind);

        pwr_sums = sum(this_pwr);
        centroid_pwr = sum(this_pwr.*this_ind) ./ pwr_sums;

        centroid_dist = sqrt(sum((centroid_ind - centroid_pwr).^2));
        
        % compute width (in both dimensions) of each detection region
        %  (simplest: fit a box around region)
        centroid_width = max(this_ind)-min(this_ind);

        %========================
        %========================
        % build output structure
        %========================
        %========================
        out.centroid_ind(ivalid_region_counter)   = centroid_ind;
        out.centroid_pwr(ivalid_region_counter)   = centroid_pwr;
        out.centroid_dist(ivalid_region_counter)  = centroid_dist;
        out.centroid_width(ivalid_region_counter) = centroid_width;
        
%         plot(centroid_pwr(2),centroid_pwr(1),'ro')
%         plot(centroid_ind(2),centroid_ind(1),'r.')
    
    end
    
end

out.nDetCentroids = ivalid_region_counter;

if(nargout == 1)
    
    varargout{1} = out;
    
elseif(nargout == 2)
    varargout{1} = out.centroid_pwr;
    varargout{2} = out.centroid_width;
    
end

end