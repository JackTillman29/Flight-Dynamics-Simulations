function [u,v,w] = uv2uvw(u,v)
    w = sqrt(1 - u.^2 - v.^2);
end