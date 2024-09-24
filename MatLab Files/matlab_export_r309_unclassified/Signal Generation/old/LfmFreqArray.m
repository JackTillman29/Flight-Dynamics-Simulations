function freq = LfmFreqArray(Fc,BW,numRamps,L)

normT = linspace(0,1,L+1);
normT = normT(1:end-1);

ramp = mod(normT,1/numRamps)*numRamps - 0.5;
freq = BW*(ramp+Fc);

end



