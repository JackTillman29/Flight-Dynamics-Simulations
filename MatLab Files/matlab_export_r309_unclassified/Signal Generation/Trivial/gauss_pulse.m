function y = gauss_pulse(t,delay,tau)

y = exp(-((t - delay)/tau).^2);

end