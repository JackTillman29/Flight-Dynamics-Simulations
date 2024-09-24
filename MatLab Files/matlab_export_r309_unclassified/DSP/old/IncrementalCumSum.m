function [x,xipeaks] = IncrementalCumSum(x,accum_thresh)
accum_flag = 0;
xipeaks = [];
if(~exist('accum_thresh','var'))
    accum_thresh = 1.0e-6;
end
for k = 2 : length(x)
    ds = x(k) - x(k-1);
    if ( ( ds > accum_thresh) && accum_flag == 0)
        accum_flag = 1;
    end
    if( ( ds < accum_thresh) && accum_flag == 1)
        accum_flag = 0;
        xipeaks = [xipeaks k-1];
        x(k:end) = x(k:end) - x(k);
    end
end
end