function [z] = s2z(s,Ts,method,Fwarp)
% This function will convert locations in the s-plane to locations in the z-plane.
%
% Once the continuous time polynomial function
% is of the form: k * (s-p1)(s-p2)(s-p3) .... (s-pn)
% 
if(~exist('method','var'))
    method = 'tustin';
    disp('Defaulting to method=tustin');
    disp('Other options: prewarp (Hz), polezero');
end

switch(lower(method))
    case 'tustin'
        Q = 2.0 / Ts;
        z = (Q+s)./(Q-s);
    case 'prewarp'
        w = 2.0 * pi * Fwarp;
        Q = w./ (tan(0.5 .* w .* Ts));
        z = (Q+s)./(Q-s);
    case 'polezero'
        z = exp(s.*Ts);
    otherwise
        error(['s2z: Unrecognized option: ' method]);
end

end