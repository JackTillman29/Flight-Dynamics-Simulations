x = ginput;
fv = 1e-3 * (sig3_fft_freq-Fif);
nb = size(x,1) - 1;
csum  = zeros(nb,1);
csumf = csum;
for k = 1 : ( size(x,1) - 1 )
    kstart = find(fv > x(k  ,1),1,'first');
    kstop  = find(fv > x(k+1,1),1,'first');
    tempsum = cumsum( sig3_psd( kstart:kstop ) ) .* df;
    csum(k) = tempsum(end);
    jmax = max(sig3_psd( kstart:kstop ));
    jmax = find(sig3_psd( kstart:kstop ) == jmax);
    csumf(k) = fv(kstart-1+jmax);
end