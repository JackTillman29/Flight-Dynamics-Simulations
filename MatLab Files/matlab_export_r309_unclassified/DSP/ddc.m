function y = ddc(u,f,Ts)
% perform direct digital down conversion
k = 0:length(u)-1;
rad_per_sample = 2*pi*f*Ts;
y = u.*exp(1i*rad_per_sample*k);

end