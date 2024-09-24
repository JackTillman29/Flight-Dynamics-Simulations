function cb_plot_single(varargin)
    hi = varargin{1};
    getpt = get(get(hi,'Parent'),'CurrentPoint');
    rownum = round(getpt(1,2));
    imdat = get(hi,'CData');
    if(isempty(get(hi,'UserData')))
        ylimits = get(gca,'CLim');
        figure;
        hp = plot(1:size(imdat,2));
%         title(['row:',num2str(rownum)])
        ylim(ylimits);
        set(hi,'UserData',hp);
    end
    
    % hp is valid now
    try
        hp = get(hi,'UserData');
        set(hp,'YData',imdat(rownum,:));
%         tmpax = get(hp,'Parent');
%         title(get(get(tmpax,'Title'),'String'),['row:',num2str(rownum)])
        
    catch
        set(hi,'UserData',[]);
    end
    
end