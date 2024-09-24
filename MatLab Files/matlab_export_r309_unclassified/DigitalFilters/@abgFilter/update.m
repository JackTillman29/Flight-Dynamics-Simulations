function f = update( f, err )
%UPDATE update ABG filter based on an externally measured error
f.xbar = f.xhat - f.gain * err;
f.xhat = f.Phi * f.xbar;
end

