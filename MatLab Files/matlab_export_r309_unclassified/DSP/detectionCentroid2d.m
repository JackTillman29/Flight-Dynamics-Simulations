function varargout = detectionCentroid2d(dets,detpwr)
% Computes the centroids of connected detection regions, which can be
% weighted by the power in each detection cell.
%
% Syntax:
% Inputs:
%   dets:              MxN binary data (common output from CFAR)
%   detpwr [optional]: MxN power data used to weight the centroid of
%                      cells in binary data "dets"
% Outputs:
%   out = detectionCentroid2d(dets,detpwr)
%       out: output structure containing the fields:
%           out.centroid_ind
%           out.centroid_pwr
%           out.centroid_dist
%           out.centroid_width
%   [region_centers, region_widths] = detectionCentroid2d(dets,detpwr)
%           region_centers: centers of connected regions (weighted by "detpwr")
%
%  Von Neumann neighborhood:
%       computes the centroids using a "Von Neumann neighborhood" around a
%       cell-under-test:
%               -   V   -
%               V  CUT  V
%               -   V   -
%       Notice the corners are not included in the test. The neighborhood
%       that does include all 8 neighbors is called a "Moore neighborhood".
%               M   M   M
%               M  CUT  M
%               M   M   M
%       The Hoshen-Kopelman Algorithm implemented here uses the Von
%       Neumann neighborhood, so a map like the following will produce 5
%       distinct detection centroids, not a single detection:
%           X  0  0  0  0  0        1  0  0  0  0  0
%           0  X  0  0  0  0        0  2  0  0  0  0
%           0  0  X  0  0  0  ==>   0  0  3  0  0  0  
%           0  0  0  X  0  0        0  0  0  4  0  0
%           0  0  0  0  X  0        0  0  0  0  5  0
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
N = size(dets);
for ix = 1:N(1)
    for iy = 1:N(2)
        if(dets(ix,iy) == 1)
            
            if(ix == 1)
                left = 0;
            else
                left = dets(ix-1,iy);
            end
            if(iy == 1)
                above = 0;
            else
                above = dets(ix,iy-1);
            end
            
            if( (left == 0)  &&  (above == 0) )
                % if not left and not above
                %   THEN NEW CLUSTER FOUND
                ilabel = ilabel + 1;
                label(ix,iy) = ilabel;
            elseif( (left == 1)  &&  (above == 0) )
                label(ix,iy) = label(ix-1,iy);
            elseif( (left == 0)  &&  (above == 1) )
                label(ix,iy) = label(ix,iy-1);
            else
                icombined = icombined + 1;
                % if left and above are already labeled, then the CUT
                % connects these two regions, so need to make sure they are
                % all labeled the same (prioritizing label assigned to the
                % LEFT cell)
                label_left = label(ix-1,iy);
                label_above = label(ix,iy-1);
                label(ix,iy) = label_left;
                % find all labels with "label_above" and set to
                % "label_left"
                tmpidx = find(label == label_above);
                label(tmpidx) = label_left;
                
%                 % just for book-keeping (not a component of
%                 % algorithm)
%                 removed_labels(icombined) = label_above;
            end
            
        else
            
        end
    end
end

maxRegionNumber = max(label(:));

% compute centroids of each region
%    to do this, compute the mean of the region's pixel locations
ix = 1:N(1);
iy = 1:N(2);
[ixg,iyg] = meshgrid(ix,iy);

ivalid_region_counter = 0;
for iregion = 1:maxRegionNumber
    
    % binary matrix with 1 where only this region exists
    % (all other regions are zeroed out)
    this_region = 1*(label == iregion);
    
    this_ind = find(label == iregion);
    
    % WORKS WITHOUT THIS IF(~ISEMPTY(this_ind))
    if(~isempty(this_ind))
        
        ivalid_region_counter = ivalid_region_counter + 1;
    
        this_ixg = ixg(this_ind);
        this_iyg = iyg(this_ind);
        this_pwr = detpwr(this_ind);
        
        % center-of-mass analogy:
        %    summation(i=1,N, (xi-CM)*mi) = 0
        %      solve for CM ("Center-of-Mass"):
        %           (x1*m1 - CM*m1) + ... + (xN*mN - CM*mN) = 0
        %           sum(i=1,N, xi*mi) = CM * sum(i=1,N, mi)
        %           CM = sum(i=1,N, xi*mi) / sum(i=1,N, mi)
        
        % if mi = 1 for all i, then CoM equation reduces to simply
        % computing the mean of the indices
        centroid_ind = [mean(this_ixg) mean(this_iyg)];

        pwr_sums = sum(this_pwr);
        centroid_pwr = [sum(this_pwr.*this_ixg) sum(this_pwr.*this_iyg)] ./ pwr_sums;

        centroid_dist = sqrt(sum((centroid_ind - centroid_pwr).^2));
        
        % compute width (in both dimensions) of each detection region
        %  (simplest: fit a box around region)
        centroid_width = [max(this_ixg)-min(this_ixg) ...
                          max(this_iyg)-min(this_iyg)];

        %========================
        %========================
        % build output structure
        %========================
        %========================
        out.centroid_ind(ivalid_region_counter,1:2)   = centroid_ind;
        out.centroid_pwr(ivalid_region_counter,1:2)   = centroid_pwr;
        out.centroid_dist(ivalid_region_counter)      = centroid_dist;
        out.centroid_width(ivalid_region_counter,1:2) = centroid_width;
        
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