function [out] = ADC_Noise(N,TotalPower)
  out = sqrt(TotalPower).*exp(1i*2*pi*rand(N,1));    
end