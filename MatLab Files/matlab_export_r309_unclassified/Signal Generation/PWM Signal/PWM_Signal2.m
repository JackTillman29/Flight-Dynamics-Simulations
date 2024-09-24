function gate = PWM_Signal2(pw,pri,dt,N,N_start,pri_mod_fcn,varargin)

t = ((0:N-1)-(N_start-1))*dt;
ramp = mod(t,pri);
gate = 1*(ramp <= pw);
gate(t < 0) = 0;

end
