function [az,el] = uvw2azel(u,v,w)
    % [az,el] = uvw2azel(u,v,w)
    % Array Face (u,v,w) definition: 
    % w = normal
    % v = up
    % u = v x w (from along normal, looking left).
    el = asin(v);
    az = atan2(-u,w);
end