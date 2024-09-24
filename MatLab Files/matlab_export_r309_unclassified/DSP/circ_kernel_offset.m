function K = circ_kernel_offset(varargin)
% Usage:
%     K = circ_kernel_offset(N, angle_deg)
%           N = [n1 n2];  OR N is scalar
%     K = circ_kernel_offset(M, N, angle_deg)
%     K = circ_kernel_offset(N, angle_deg, dist_from_ctr, circ_radius)
%           N = [n1 n2];  OR N is scalar
%     K = circ_kernel_offset(M, N, angle_deg, dist_from_ctr, circ_radius)

dist_from_center = 0.5;
circ_radius      = 0.5;
if(nargin == 2)
    % K = circ_kernel_offset(N, angle_deg)
    if(length(varargin{1}) == 1)
        M = varargin{1};
        N = varargin{1};
    elseif(length(varargin{1}) == 2)
        tmp = varargin{1};
        M = tmp(1);
        N = tmp(2);
    else
        error('first input variable must be length-1 or length-2')
    end
    angle_deg = varargin{2};
elseif(nargin == 3)
    % K = circ_kernel_offset(M, N, angle_deg)
    M = varargin{1};
    N = varargin{2};
    angle_deg = varargin{3};
elseif(nargin == 4)
    % K = circ_kernel_offset(N, angle_deg, dist_from_ctr, circ_radius)
    %       N = [n1 n2];
    if(length(varargin{1}) == 1)
        M = varargin{1};
        N = varargin{1};
    elseif(length(varargin{1}) == 2)
        tmp = varargin{1};
        M = tmp(1);
        N = tmp(2);
    else
        error('first input variable must be length-1 or length-2')
    end
    angle_deg = varargin{2};
    dist_from_center = varargin{3};
    circ_radius = varargin{4};
elseif(nargin == 5)
    % K = circ_kernel_offset(M, N, angle_deg, dist_from_ctr, circ_radius)
    M = varargin{1};
    N = varargin{2};
    angle_deg = varargin{3};
    dist_from_center = varargin{4};
    circ_radius = varargin{5};
else
    error('wrong number of inputs')
end
% % % 
% % % x = linspace(-0.5,0.5,N);
% % % y = linspace(-0.5,0.5,M);
% % % 
% % % [X,Y] = meshgrid(x,y);
% % % 
% % % K = zeros(M,N);
% % % 
% % % K(sqrt(X.^2 + Y.^2)<=0.5) = 1;

cx = dist_from_center*cosd(angle_deg);
cy = dist_from_center*sind(angle_deg);

x = linspace(-dist_from_center-circ_radius, dist_from_center+circ_radius,N);
y = linspace(-dist_from_center-circ_radius, dist_from_center+circ_radius,M);

[X,Y] = meshgrid(x,y);

K = zeros(M,N);

K(sqrt((X-cx).^2 + (Y-cy).^2) <= circ_radius) = 1;


end