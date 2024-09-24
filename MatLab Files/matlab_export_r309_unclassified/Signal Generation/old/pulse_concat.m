function out = pulse_concat(input_pri,n,phaseShift)
    if ( ~exist('phaseShift') ) 
        phaseShift = zeros(1,n);
    else
        % test length of phase shift
        if(length(phaseShift) ~= n)
            error(['pulse_concat needs ' num2str(n) ' phase shift values. You only provided ' num2str(length(phaseShift)) '.']);
        end
    end

    % build the command
    s = 'out = [ ';
    for k = 1:n
        s = [s 'exp(1j*phaseShift(' num2str(k) '))*input_pri '];
    end
    s = [s '];'];
    eval(s);
end