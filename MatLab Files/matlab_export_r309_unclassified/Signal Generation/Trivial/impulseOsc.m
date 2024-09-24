function x = impulseOsc(f,T,dt)

x = clockOsc(f,T,0.5,dt);
dx = diff([0 x]);

x = zeros(size(x));
x(dx > 0) = 1;

end