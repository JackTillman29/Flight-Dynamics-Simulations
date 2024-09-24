function v = PWM_Signal(pw,pri,dt,N,N_start,pri_mod_fcn,varargin)
% function v = PWM_Signal( ...
%        pw, ...    % pulsewidth (sec)
%        pri, ...   % pri (sec)
%        dt, ...    % sample time (sec)
%        N, ...     % # of samples to make signal
%        N_start, ... % [opt] sample to begin first rise
%        pri_mod_fcn) % [opt] pri modulation function handle
%        varargin     % input modulation function argumens
%
%                     % NOTE: modulation function expects input as
%                     %       fraction of total length.
%                     %       Expects output as % change from input

if(~exist('N_start'))
    N_start = 1;
end

if(~exist('pri_mod_fcn'))
    pri_mod_fcn = [];
end

if(isempty(pri_mod_fcn))
    
    nOn  = round(pw/dt);
    nOff = round(pri/dt) - nOn;
    nSingle = nOn + nOff;
    
    x = zeros(1,nSingle);
    x(1:nOn) = 1;
    
    nLoops = ceil(N/nSingle);
    
    concat_start    = 1;
    concat_stop     = concat_start + nSingle - 1;
    
    v = zeros(1,nLoops*nSingle);
    
    for k = 1 : nLoops
        v(concat_start:concat_stop) = x;
        concat_start = concat_stop + 1;
        concat_stop  = concat_start + nSingle - 1;
    end
    
    v = v(1:N);
    v = [zeros(1,N_start-1) v(1:(N-N_start+1))];
    
else
    v = zeros(1,N);
    cFrac = N_start/N;
    idx1 = N_start;
    idx2 = N_start;
    done = 0;
    while(~done)
        [pri_fac,pw_fac]=pri_mod_fcn(idx2/N,varargin{:});
        % get first pulsewidth -- untouched
        pw_prime = pw * pw_fac;
        
        
        idx2 = idx1 - 1 + round(pw_prime / dt);
        if(idx2 >= N)
            idx2 = N;
            done = 1;
        end
        
        % set signal high
        v(idx1:idx2) = 1;
        
        % determine next rising edge using function
        pri_prime = pri * pri_fac;
        
        idx1 = idx2 + max(0,round((pri_prime-pw_prime) / dt)) + 1;

    end
    
end
end