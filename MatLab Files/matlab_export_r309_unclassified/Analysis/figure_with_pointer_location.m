function figure_with_pointer_location()

hfig = figure;

set(hfig,'WindowButtonMotionFcn',@mouseMove)

% hax = axes;
data.htb = uicontrol('Style','Text', 'Units','Normalized', ...
    'Position',[0 0 1.0 0.05], ...
    'HorizontalAlignment','right',...
    'String','this is a text box');

guidata(hfig,data)

    function mouseMove(obj,eventdata)
        data = guidata(obj);
        fmtstr = '%10.5f';
        pt = get(gca,'CurrentPoint');
        set(data.htb,'String',['(x,y) = (',num2str(pt(1,1),fmtstr),', ',num2str(pt(1,2),fmtstr),')'])
            
%         end
        
        guidata(obj,data);
    end


end