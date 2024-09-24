function figprint(fnum,inOutputName,inOutputRes,inAxesSize,inLabelSize,inTitleSize,inLineSize,inGridLineStyle)
% default: figprint(2,'LFM_Advance1p2Mhz.png',300,8,12,12,1.5,'-')
    % determine the subplots
    subplots = findobj('Type','axes','Parent',fnum);
    xlabels = [];
    ylabels = [];
    zlabels = [];
    titles  = [];
    lines   = [];
    % loop through all subplots
    for ksubplot=1:length(subplots)
        thisaxis = subplots(ksubplot);
        xlabels = [xlabels get(thisaxis,'XLabel')];
        ylabels = [ylabels get(thisaxis,'YLabel')];
        zlabels = [zlabels get(thisaxis,'ZLabel')];
        titles  = [titles  get(thisaxis,'Title')];
        lines   = [lines findobj('Type','line','Parent',thisaxis)'];
        set(thisaxis,'GridLineStyle',inGridLineStyle);
        box on;
    end % loop through all subplots
    set(xlabels,'FontSize',inLabelSize);
    set(ylabels,'FontSize',inLabelSize);
    set(titles,'FontSize',inTitleSize);
    set(subplots,'FontSize',inAxesSize);
    set(lines,'LineWidth',inLineSize);
    
    print('-dpng',['-f' num2str(fnum)],['-r' num2str(inOutputRes)],inOutputName);
end