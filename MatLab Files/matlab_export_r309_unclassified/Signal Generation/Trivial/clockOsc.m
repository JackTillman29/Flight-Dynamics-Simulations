function x = clockOsc(f,T,dutyCycle,dt)
% x = clockOsc(f,T,dutyCycle,dt)
%
% tip: to get an impulse oscillator (one high sample per period), set
% dutyCycle equal to dt:
%     impulseOsc = clockOsc(f,T,dt,dt)
% 

if(isinteger(T))
    N = double(T);
else
    N = floor(T/dt);
end

t = [0:N-1]*dt;
period = 1/f;

tmod = mod(t,period)./period;

x = zeros(1,N);
x(find(tmod < dutyCycle)) = 1;

end