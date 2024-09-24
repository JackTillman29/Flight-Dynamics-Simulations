function wo_t_f_plot(wo,fig,descr)
% takes an existing waveobject plot and re-formats to custom time/freq
% domain plot
    ah=findobj('Type','axes','Parent',fig);
    linkaxes(ah,'off');
    plot(ah(1), 1e6*(0:length(wo.tdSignal)-1) ./ wo.Fs, ...
        1e3*real(wo.tdSignal));
    grid on;
    xlabel('Time (\mus)');
    ylabel('Amplitude (mV)');
    title('Time Domain Signal');
    title(ah(2), ...
        {descr,['Spectrum Power (\Sigma=' num2str(10*log10(wo.cumpow),'%8.1f') 'dBW)']});
    set(ah(1),'Tag','testpoint1');
    set(ah(2),'Tag','testpoint2');
end