function p = patchline(xs,ys,varargin)
% syntax
%       patchline(xs,ys)
%       patchline(xs,ys,zs)
%       patchline(xs,ys,'PropName',propval)
%       patchline(xs,ys,zs,'PropName',propval)

[zs,PVs] = parseInputs(varargin{:});
if(mod(numel(PVs),2) ~= 0)
    % odd number of inputs
    error('patchline: Param-Value must be entered in valid pairs')
end

if(isempty(zs))
    p = patch([xs(:);NaN],[ys(:);NaN],'k');
else
    p = patch([xs(:);NaN],[ys(:);NaN],[zs(:);NaN],'k');
end


% apply PV pairs
for k = 1:2:numel(PVs)
    set(p,PVs{k},PVs{k+1});
end
if nargout == 0
    clear p
end

    function [zs,PVs] = parseInputs(varargin)
        if(isnumeric(varargin{1}))
            zs = varargin{1};
            PVs = varargin{2:end};
        else
            PVs = varargin;
            zs = [];
        end
    end

end