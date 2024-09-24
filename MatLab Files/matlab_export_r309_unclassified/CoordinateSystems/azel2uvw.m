function [u,v,w] = azel2uvw(azrad,elrad)
    % function [u,v,w] = azel2uvw(azrad,elrad)
    % Array Face (u,v,w) definition: 
    % w = normal
    % v = up
    % u = v x w (from along normal, looking left).
    u = -cos(elrad).*sin(azrad);
    v = sin(elrad);
    w = cos(elrad).*cos(azrad);
end