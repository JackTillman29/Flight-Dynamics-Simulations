function filt = ChebyshevFilterBandPass(filcen,filbw,ripplepct,npoles,Fs)
% copy : ChebyshevFilterBandPass(filcen,filbw,ripplepct,npoles,Fs)
% function filt = ChebyshevFilterBandPass(
% filcen     Hz
% filbw      Hz
% ripplepct  0-29
% npoles     2-20
% Fs         Hz
% NOTE: No output arguments presents plot

FC_LO = (filcen + filbw/2)/Fs;
FC_HI = (filcen - filbw/2)/Fs;
PCT_RPL = ripplepct;
NPOLES = npoles;
[filt.num_lo,filt.den_lo]=ChebyshevFilter(FC_LO,0,PCT_RPL,NPOLES);
[filt.num_hi,filt.den_hi]=ChebyshevFilter(FC_HI,1,PCT_RPL,NPOLES);

if(nargout == 0)
    disp(['PCT: ' num2str([FC_LO FC_HI])]);
    Ts = 1/Fs;
    tdur = 5e-3;
    bw_plot_span = 200;
    t = 0:Ts:tdur;
    xi = (0:(length(t)-1))./length(t);
    fi = bw_plot_span*filbw.*xi-0.5*bw_plot_span*filbw;
    figure;
    uin=SimpleLfm(filcen-bw_plot_span*filbw,filcen+bw_plot_span*filbw,t(end),1/Fs);
    uoutlo = applyIIR(filt.num_lo,filt.den_lo,uin,'double');
    uouthi = applyIIR(filt.num_hi,filt.den_hi,uin,'double');
    uout   = applyIIR(filt.num_hi,filt.den_hi,uoutlo,'double');
    plot(1e-3*fi,20*log10(abs([uoutlo.' uouthi.' uout.'])));
    xlabel(['kHz rel ' num2str(1e-3*filcen) 'kHz']);
    ylabel('20Log_{10}(A)');
    ylim([-50 5]);
    %xlim(1e-3*filbw*[-4 4]);
    title('Filter Response Curve');
    grid on;
    
    try
        pp = PrepForPrint();
    catch
        addpath('z:\sawmillkd\MATLAB\Printing');
        pp = PrepForPrint();
    end
    add_print_callbacks;
    
end

end