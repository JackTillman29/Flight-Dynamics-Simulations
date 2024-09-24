function [y] = applyFilter( f, u )
%step stimulate ABG with a step input. f,tstart,tdur. tstart is the step time
y(1,:) = [0 f.xbar' f.xhat'];

for k = 2 : length(u)
    % NOTE: Error term is ESTIMATE - MEASUREMENT
    err = f.xhat(1)-u(k);
    f = update(f,err);
    y(k,:) = [err f.xbar' f.xhat'];
end

if(nargout == 0)
    figure(1);
    
    subplot(4,1,1);
    plot([y(:,2) u'],'.-');
    grid on; xlabel('Time (sec)');  ylabel('x');
    
    subplot(4,1,2);
    plot([y(:,3)],'.-');
    grid on; xlabel('Time (sec)');  ylabel('xdot');
    
    subplot(4,1,3);
    plot([y(:,4)],'.-');
    grid on; xlabel('Time (sec)');  ylabel('xdot');
    
    subplot(4,1,4);
    plot([y(:,1)],'.-');
    grid on; xlabel('Time (sec)');  ylabel('err');
    
end

end

