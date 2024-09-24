function cm = colormap_fade2black(varargin)
% colormap_fade2black(varargin)
%    varargin{1} = N        : number of colors in colormap
%    varargin{2} = fraction : fraction of colormap that is faded to black
%                  [0-1]
if(nargin == 0)
    N        = 256;
    fraction = 0.5;
    refmap = @jet;
elseif(nargin == 1)
    N        = varargin{1};
    fraction = 0.5;
    refmap = @jet;
elseif(nargin == 2)
    N        = varargin{1};
    fraction = varargin{2};
    refmap = @jet;
elseif(nargin == 3)
    N        = varargin{1};
    fraction = varargin{2};
    refmap   = varargin{3};
end

Nb = round(N*fraction);
N = N - Nb;
cm = refmap(N);
cmb = linspace(0,0.51,Nb);
cm = [0*cmb' 0*cmb' cmb';cm];

end