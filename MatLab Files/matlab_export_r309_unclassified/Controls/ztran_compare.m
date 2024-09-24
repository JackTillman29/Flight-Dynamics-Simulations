Ts = 0.001;

s = rand() + 1j * rand()

z1 = (1 + (Ts/2).*s)./(1-(Ts/2).*s)
z2 = exp(s.*Ts)