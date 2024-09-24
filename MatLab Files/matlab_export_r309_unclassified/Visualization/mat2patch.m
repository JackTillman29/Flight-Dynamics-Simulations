function [p_x, p_y, p_c, varargout] = mat2patch(mat,varargin)
% Inputs:
%   mat: MxN matrix of values
%           for example, "mat" out of imread
%           [mat,mat_cmap] = imread(img_filename);
%   alpha_map [optional]: MxN matrix of alpha values (between 0 and 1)
%
% Outputs:
%   p_x: vector of x locations (4xN)
%   p_y: vector of y locations (4xN)
%   p_c: vector of c locations (1xN)
%   varargout{1}: alpha map of each patch (1xN)
%
% % author: Jeff Hole (Booz Allen Hamilton), 2019-12

if(length(varargin) == 1)
    alpha_map = varargin{1};
else
    alpha_map = ones(size(mat));
end

% if mat is a MxNx3, it is a matrix of RGB triplets
mat_is_rgb_triplets = 0;
if(length(size(mat)) == 3)
    mat_is_rgb_triplets = 1;
end

% if(length(varargin) > 0)
%     xv = varargin{1};
%     yv = varargin{2};
%     if(length(varargin) > 2)
%         zv = varargin{3};
%     end
% else
%     xv = 1:size(mat,1);
%     yv = 1:size(mat,2);
% %     zv = 1;
% end

Lx = size(mat,1);
Ly = size(mat,2);

dx = 1;
dy = 1;

npts = prod(size(mat));
p_x  = zeros(4,npts);
p_y  = zeros(4,npts);
% p_z  = zeros(4,npts);
if(mat_is_rgb_triplets == 0)
    p_c = zeros(npts,1);
else
    p_c = zeros(npts,3);
end
p_alpha = zeros(npts,1);

ipt = 0;
for ix = 1:Lx
    for iy = 1:Ly
        ipt = ipt + 1;
        ctr = [ix iy];
        p_x(:,ipt) = [ix-dx/2 ix+dx/2 ix+dx/2 ix-dx/2].' / Lx;
        p_y(:,ipt) = [iy-dy/2 iy-dy/2 iy+dy/2 iy+dy/2].' / Ly;
%         p_z(:,ipt) = 1 * ones(4,1);
        if(mat_is_rgb_triplets == 0)
            p_c(ipt,1) = 1.0*mat(ix,iy);
        else
            p_c(ipt,1:3) = single(squeeze(mat(ix,iy,:)))/255;
        end
        
        p_alpha(ipt,1) = alpha_map(ix,iy);
    end
end

p_x = p_x - 0.5;
p_y = p_y - 0.5;

if(nargout > 3)
    varargout{1} = p_alpha;
end


end