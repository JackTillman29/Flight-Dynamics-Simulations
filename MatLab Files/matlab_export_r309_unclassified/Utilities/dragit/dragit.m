% written by Keith Sawmiller and Jeff Hole

function varargout = dragit()
    
    hfig = gcf;
    hax = gca;
    
    % get current (last) line object
    hchild = get(gca,'Children');
    
    hL  = hchild(end);
    
    hL_ref = copyobj(hL,hax);
    set(hL_ref,'Color','r','LineStyle',':');
    
    
    % set line objects as UserData of the axes object
    userDat = get(gcf,'UserData');
    userDat.dragit.hfig   = hfig;
    userDat.dragit.hax    = hax;
    userDat.dragit.hL     = hL;
    userDat.dragit.hL_ref = hL_ref;
    userDat.dragit.data = [0 0; 0 0];
    userDat.dragit.hmsg = msgbox('When you''re done, press any key...')
    set(hfig,'UserData',userDat);
    
%     set(gcf,'UserData',[0 0; 0 0]);
    set(hfig,'WindowButtonDownFcn',@cb_drag_down);
    set(hfig,'WindowButtonUpFcn',@cb_drag_up);
    set(hfig,'WindowKeyReleaseFcn',@cb_dragit_done);
    grid on;
    
    if(nargout == 1)
        varargout{1} = hL;
    end
    
    
end


function cb_dragit_done(varargin)

    userDat = get(gcbo,'UserData');
    hL = userDat.dragit.hL;
    xsave = get(hL,'XData');
    ysave = get(hL,'YData');
    xout = [xsave',ysave'];
    
    % delete reference line
    delete(userDat.dragit.hL_ref);
    
    % prompt user to save data to file (will save to working directory)
    prompt = {'Enter output text filename:'};
    dlg_title = 'Save modified data?';
    num_lines = 1;
    def = {''};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if(~isempty(answer))
        hmsg = msgbox(['Data saved to ',answer])
        dlmwrite(answer{1},xout,'\n')
    else
        hmsg = msgbox(['Did not save modified data!'])
    end
    
    % remove callbacks for dragit from the figure
    set(userDat.dragit.hfig,'WindowButtonDownFcn','');
    set(userDat.dragit.hfig,'WindowButtonUpFcn','');
    set(userDat.dragit.hfig,'WindowKeyPressFcn','');
    
    close(userDat.dragit.hmsg)

end


function cb_drag_down(varargin)
    current_pt = get(gca,'CurrentPoint');
    userDat = get(gcbo,'UserData');
    Q = userDat.dragit.data;
    Q(1,:) = current_pt(1,1:2);
    userDat.dragit.data = Q;
    set(gcbo,'UserData',userDat);
end


function cb_drag_up(varargin)
    current_pt = get(gca,'CurrentPoint');
    userDat = get(gcbo,'UserData');
    Q = userDat.dragit.data;
    Q(2,:) = current_pt(2,1:2);
    userDat.dragit.data = Q;
    set(gcbo,'UserData',userDat);
    
%     hL = evalin('caller','hL');
    hL = userDat.dragit.hL;
    
    % Have points, now find actual points
    xd = get(hL,'XData');
    yd = get(hL,'YData');
    
    rd_1 = sqrt( (xd-Q(1,1)).^2 + (yd-Q(1,2)).^2 );
    % seems to work better using x distance only
%     rd_1 = sqrt((xd-Q(1,1)).^2);
    
    [m1,idx1] = min(rd_1);
    
    xd(idx1) = Q(2,1);
    yd(idx1) = Q(2,2);
    
    set(hL,'XData',xd,'YData',yd);
    assignin('base','u',xd);
    assignin('base','y',yd);
end