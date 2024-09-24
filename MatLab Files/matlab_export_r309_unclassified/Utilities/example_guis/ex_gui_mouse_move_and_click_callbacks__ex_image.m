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
            
            if(~isfield(data,'ipt'))
                data.ipt = 0;
                data.pts = 0;
            end
            data.ipt = data.ipt + 1;
            data.pts(data.ipt,1:2) = pt(1,1:2);
            
            hax1 = subplot(3,1,1:2);
            plot(data.pts(:,1),data.pts(:,2),'.')
            xlim([-10 10])
            ylim([-10 10])
            title(gca,...
                {['points collected: ',num2str(data.ipt)],...
                ['(x,y) = (',num2str(pt(1,1)),', ',num2str(pt(1,2)),')']})
            
            hax2 = subplot(3,1,3);
            plot(data.pts(:,1))
            hold on;
            plot(data.pts(:,2))
            hold off;
            
            axes(hax1)
            
            % save each point
            
        else
            if(isfield(data,'ipt'))
                data.ipt = 0;
                data.pts = 0;
            end
        end
        
        guidata(obj,data);
    end


end