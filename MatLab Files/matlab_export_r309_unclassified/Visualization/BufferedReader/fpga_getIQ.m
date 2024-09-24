function y = fpga_getIQ(u,varargin)
% Purpose of this file is to package FPGA DDR3 block
% data writes in the SD card design
    % varargin{1} = bias
    % varargin{2} = scale factor
    if( isempty(varargin{:}) )
        bias = 0.0;
        scale_factor = 1.0;
    else
        bias = varargin{1}{1};
        scale_factor = varargin{1}{2};
    end
    y =      scale_factor * single(u(1:2:end)) + bias + ...
       1i * (scale_factor * single(u(2:2:end)) + bias);
end