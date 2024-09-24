function K = circ_kernel(varargin)

if(nargin == 1)
    if(length(varargin{1}) == 1)
        M = varargin{1};
        N = varargin{1};
    elseif(length(varargin{1}) == 2)
        tmp = varargin{1};
        M = tmp(1);
        N = tmp(2);
    else
        error('input variable must be length-1 or length-2')
    end
elseif(nargin == 2)
    M = varargin{1};
    N = varargin{2};
else
    error('wrong number of inputs')
end

x = linspace(-0.5,0.5,N);
y = linspace(-0.5,0.5,M);

[X,Y] = meshgrid(x,y);

K = zeros(M,N);

K(sqrt(X.^2 + Y.^2)<=0.5) = 1;


end