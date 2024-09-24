function ProgressBar(i,N,varargin)
% TimeMe(i,N,varargin)
%   i : current loop index
%   N : number of loop iterations
%   varargin{1} : print progress message every "varargin{1}" loops
%   varargin{2} : <string> message string to print before percentage.


% defaults
skip = 1;
message = '';

% 3rd arg: number of 
% 4th arg: message string
if(nargin > 2)
    if(nargin == 3)
        skip = varargin{1};
    elseif(nargin == 4)
        skip = varargin{1};
        message = varargin{2};
    end
end

if(mod(i,skip) == 0)
clc
disp([message '.....' num2str(floor(i*100/N)) '%'])
end

end