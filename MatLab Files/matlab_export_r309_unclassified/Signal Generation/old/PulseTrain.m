function out = PulseTrain(prf,pw,N,dt,varargin)
% Help: PulseTrain(prf,pw,N,dt)
%       prf = pulse repitition frequency (Hz)
%       pw  = pulsewidth (sec)
%       N   = total length of output vector (# elements)
%       dt  = sample time (sec)
%       {1} -> 'length' defines N as stated above
%              '#prix' defines N as number of pris
    pri = 1.0 / prf;
    if(nargin > 4)
        switch(lower(varargin{1}))
            case 'length'
               % do nothing
            case '#pris'
                N = floor((N*pri)/(dt));
            otherwise
                error('unknown option in PulseTrain');
        end
    end
    onTime  = round(pw / dt);
    offTime = round( (pri - pw) / dt );
    totTime = onTime + offTime;
    iCount = ceil( N/totTime );
    out = dt * (0:(N-1));
    for k = 1 : iCount
        jstart = (k-1) * totTime + 1;
        jstop  = jstart + onTime - 1;
        kstart = jstop + 1;
        kstop  = kstart + offTime - 1;

        out(jstart:jstop) = 1;
        out(kstart:kstop) = 0;
    end
    out = out(1:N);
end