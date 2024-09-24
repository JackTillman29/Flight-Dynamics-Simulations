function ex_gui_mouse_follow()

close all;
hfig = figure;

set(hfig,'WindowButtonMotionFcn',@mouseMove)
set(hfig,'WindowButtonDownFcn',@mouseClickDown)
set(hfig,'WindowButtonUpFcn',@mouseClickUp)

hax = axes;

%     function mouseMove(obj,eventdata)
% 
%         c = get(gca,'CurrentPoint');
% 
%         eventdata
% 
%         title(gca,['(x,y) = (',num2str(c(1,1)),', ',num2str(c(1,2)),')'])
% 
%     end

data.mouseClickStatus = 0;

guidata(hfig,data)


    function mouseClickDown(obj,eventdata,varargin)
        data = guidata(obj);
        %eventdata
        data.mouseClickStatus = 1;
        guidata(obj,data);
        
    end

    function mouseClickUp(obj,eventdata)
        data = guidata(obj);
        %eventdata
        data.mouseClickStatus = 0;
        guidata(obj,data);
    end

    function mouseMove(obj,eventdata)
        data = guidata(obj);
        %eventdata
        
        clc
        if(data.mouseClickStatus == 1)
            disp('mouse button down')
            % PUT ACTIONS YOU WANT TO TAKE HERE
            % like getting the point on the axes
            pt = get(gca,'CurrentPoint');
            
            
            plot(pt(1,1),pt(1,2),'.')
            xlim([-10 10])
            ylim([-10 10])
            title(gca,['(x,y) = (',num2str(pt(1,1)),', ',num2str(pt(1,2)),')'])
            
            
        end
        
        guidata(obj,data);
    end


end