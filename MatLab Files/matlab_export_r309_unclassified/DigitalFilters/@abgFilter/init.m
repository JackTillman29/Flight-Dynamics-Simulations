function f = init( f,xi )
%INIT initialize xbar
%   xhat is computed based on sample time
f.xbar = xi;
f.xhat = f.Phi * f.xbar;


end

