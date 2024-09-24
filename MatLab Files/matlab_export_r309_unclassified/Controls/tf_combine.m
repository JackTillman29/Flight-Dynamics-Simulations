function [out] = tf_combine(tf)
    tf_zeros = [];
    tf_poles = [];
    tf_gain  = 1;

    for i = 1:length(tf)
        [z,p,k] = zpk(tf(i).num,tf(i).den);
        tf_zeros = [ tf_zeros z ];
        tf_poles = [ tf_poles p' ];
        tf_gain  = tf_gain * k;
    end

    if(isempty(tf_zeros))
        tf_zeros = [];
    end
    out.num = tf_gain .* poly(tf_zeros);
    out.den =            poly(tf_poles);
end