function cb_move(varargin)
% src figure and dest figure need to have the same axes (subplot) layout to
% work as intended.

% get figure handle
hsrcfig = gcf;

% get axes handles
hsrcaxes = findobj('Parent',hsrcfig,'Type','axes');

% popup window with figure number input to
% get handle to final figure
hdestfig = inputdlg('Enter figure to move data to...');
hdestfig = str2num(hdestfig{1});
hdestfig = figure(hdestfig);
% get axes handles on final figure
hdestaxes = findobj('Parent',hdestfig,'Type','axes');

% TODO: if multiple axes are in the window, issue another dialog box with
% subplot numbers:
if(length(hdestaxes) > 1)
    for k = 1:length(hdestaxes)
        haxpos = get(hdestaxes(k),'Position')
        hanno(k) = annotation(hdestfig,'textbox',haxpos, ...
            'String',num2str(k),'FontSize',32,'BackgroundColor','w','FaceAlpha',0.5,'EdgeColor','none',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    idx_hdestaxes = inputdlg('Enter subplot...');
    idx_hdestaxes = str2num(idx_hdestaxes{1});
    hdestaxes = hdestaxes(idx_hdestaxes);
    delete(hanno);
end

% count the number of lines on destination axes
hp = findobj('Parent',hdestaxes,'Type','line');

a = 1;
inputcolor = 'b';
outputcolor = 'r';
for i = 1:length(hdestaxes)
    hp = findobj('Parent',hsrcaxes(i),'Type','line');
    hp_new = copyobj(hp,hdestaxes(i))
    set(hp_new,'Color','r')
    
%     % get CPSD power in title of src figure
%     pwrstr = get(get(hsrcaxes(i),'Title'),'String');
%     i1 = strfind(pwrstr,'[');
%     i2 = strfind(pwrstr,']');
%     if(~isempty(i1))
%         pwrstr_src = pwrstr(i1:i2);
%         is_hax = 1;
%     else
%         pwrstr_src = [];
%         is_hax = 0;
%     end
%     
%     % get CPSD power in title of dest figure
%     pwrstr = get(get(hdestaxes(i),'Title'),'String');
%     i1 = strfind(pwrstr,'[');
%     i2 = strfind(pwrstr,']');
%     if(~isempty(i1))
%         pwrstr_dest = pwrstr(i1:i2);
%     else
%         pwrstr_dest = [];
%     end
%     
%     % add a legend to the dest figure
%     if(is_hax)
%         legend(hdestaxes(i),{pwrstr_dest,pwrstr_src})
%     end
end

end
