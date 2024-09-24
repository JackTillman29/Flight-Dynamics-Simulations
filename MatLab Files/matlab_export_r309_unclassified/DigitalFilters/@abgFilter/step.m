function [y,t] = step( f, tstart, tdur )
%step stimulate ABG with a step input. f,tstart,tdur. tstart is the step time
dt = f.Phi(1,2);
t = 0:dt:tdur;
y = zeros(length(t),7);
u = t > tstart;

init(f,[0 0 0]');

y(1,:) = [0 f.xbar' f.xhat'];

for k = 2 : length(t)
    % NOTE: Error term is ESTIMATE - MEASUREMENT
    err = f.xhat(1)-u(k);
    f = update(f,err);
    y(k,:) = [err f.xbar' f.xhat'];
end

if(nargout == 0)
    figure;
    
    subplot(4,1,1);
    plot(t,[y(:,2) u'],'.-');
    grid on; xlabel('Time (sec)');  ylabel('x');
    
    subplot(4,1,2);
    plot(t,[y(:,3)],'.-');
    grid on; xlabel('Time (sec)');  ylabel('xdot');
    
    subplot(4,1,3);
    plot(t,[y(:,4)],'.-');
    grid on; xlabel('Time (sec)');  ylabel('xdot');
    
    subplot(4,1,4);
    plot(t,[y(:,1)],'.-');
    grid on; xlabel('Time (sec)');  ylabel('err');
    
end

end

