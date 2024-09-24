function [hax] = MakeSubplots(xData,yDataCell,yDataLabelCell,figureSize,...
                           xLabelString,titleString)

% xData       : an array of x data
% yDataCell   : a cell array of the different data to plot
% yDataLabels : a cell array of labels for each corresponding yDataCell
%               element
% figureSize  : 'small', 'medium', 'large', 'wide'



CreateFigure(figureSize)
add_print_callbacks
hzoom=zoom;
hzoom.ActionPostCallback = @callback_update_xaxis_with_timestamps
hdatacursor = datacursormode(gcf);
hdatacursor.UpdateFcn = @DataCursor_FormatTimeStamp;

numPlots = length(yDataCell);
plotCount = 0;

for i = 1:numPlots
    plotCount = plotCount + 1;
    hax(plotCount) = subplot(numPlots,1,plotCount)
    plot(xData,yDataCell{i})
    ylabel(yDataLabelCell{i})
    
    
    update_xaxis_with_timestamp(gca)
    
    
    if(i == 1)
        title(titleString)
    end
    
    if(i == numPlots)
        xlabel(xLabelString)
    end
    
end

linkaxes(hax,'x')

end