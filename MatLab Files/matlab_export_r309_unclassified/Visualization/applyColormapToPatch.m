function varargout = applyColormapToPatch(cm,data,data_limits)
%   applyColormapToPatch(cm, data, data_limits)
%       cm: a Nx3 (RGB) matrix
%       data:
%           - patch handle
%           - patch data
%       data_limits [optional]
%           - if not present, use min/max of "data" if  it is not a handle
%
% author: Jeff Hole (Booz Allen Hamilton), 2019-12

if(ishandle(data))
    hpatch = data;
    data = get(hpatch,'cdata');
    % if cdata is an Nx1x3, it is RGB triplets array
    if(size(data,3) == 3)
        ud = get(hpatch,'UserData');
        data = ud.orig_patch_data;
    else
        % add the original data to the 
        ud = get(hpatch,'UserData');
        ud.orig_patch_data = data;
        set(hpatch,'UserData',ud)
    end
end

if( ~exist('data_limits') )
    data_limits = [min(data) max(data)];
end
ncmPts = size(cm,1);
pcm_ind = zeros(ncmPts,1);
% clip data to min/max before computing indices (otherwise, interp1 will
% return NaNs for points outside the range)
data(data > data_limits(2)) = data_limits(2);
data(data < data_limits(1)) = data_limits(1);
pcm_ind = round(interp1(data_limits,[1 ncmPts],data,'linear'));
pcm = reshape(cm(pcm_ind,:),[length(pcm_ind),1,3]);

if(exist('hpatch'))
    set(hpatch,'CData',pcm);
end

if(nargout > 0)
    varargout{1} = pcm;
end

end