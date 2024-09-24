function scroll_fcn(varargin)
ch=get(varargin{1},'Children');
%varargin{2}.VerticalScrollAmount
    for k = 2 : 6
        cp=get(ch(k),'CurrentPoint');
        if(abs(cp(1,2)) < 1.0)
            active_channel = 7-k;
            
            if(varargin{2}.VerticalScrollCount == -1)
                %disp('yes');
                evalin('base',sprintf('p(%d) = p(%d) + rstep;',active_channel,active_channel));
            else
                %disp('no');
                evalin('base',sprintf('p(%d) = p(%d) - rstep;',active_channel,active_channel));
            end
            
            ac = active_channel;
            
            hp = evalin('base','hp;');
            ht = evalin('base','ht;');
            p = evalin('base','p;');
            
            ydata = evalin('base',sprintf('yfcn(f(%d),p(%d),t);',ac,ac));
            set(hp(ac),'YData',ydata);
            set(ht(ac),'String',sprintf('%1.2f',p(ac)));
            
            
            ydata_full = evalin('base',sprintf('yfcnc(f,p,t);'));
            set(hp(6),'YData',ydata_full);
            set(ht(6),'String',sprintf('Peak: %1.2f',max(ydata_full)));
        end
    end
end

%hp(k) = plot(yfcn(f(k),p(k),t));ht(k) = title(sprintf('%1.2f',p(k)));