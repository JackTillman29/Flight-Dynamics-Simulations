function cb_move_waveObj(varargin)
% src figure and dest figure need to have the same axes (subplot) layout to
% work as intended.

% get figure handle
hsrcfig = gcf;
userData = get(hsrcfig,'UserData');
wavObjSrc = userData.waveObjs;

% get axes handles
hsrcaxes = findobj('Parent',hsrcfig,'Type','axes');

% popup window with figure number input to
% get handle to final figure
hdestfig = inputdlg('Enter figure to move data to...');
hdestfig = str2num(hdestfig{1});
% get axes handles on final figure
hdestaxes = findobj('Parent',hdestfig,'Type','axes');

% get destination wave object(s)
userData = get(hdestfig,'UserData');
wavObjDest = userData.waveObjs;

plot([wavObjSrc wavObjDest])

% % a = 1;
% % inputcolor = 'b';
% % outputcolor = 'r';
% % for i = 1:length(hdestaxes)
% %     hp = findobj('Parent',hsrcaxes(i),'Type','line');
% %     set(hp,'Color','r')
% %     set(hp,'Parent',hdestaxes(i));
% %     
% %     % get CPSD power in title of src figure
% %     pwrstr = get(get(hsrcaxes(i),'Title'),'String');
% %     i1 = strfind(pwrstr,'[');
% %     i2 = strfind(pwrstr,']');
% %     if(~isempty(i1))
% %         pwrstr_src = pwrstr(i1:i2);
% %         is_hax = 1;
% %     else
% %         pwrstr_src = [];
% %         is_hax = 0;
% %     end
% %     
% %     % get CPSD power in title of dest figure
% %     pwrstr = get(get(hdestaxes(i),'Title'),'String');
% %     i1 = strfind(pwrstr,'[');
% %     i2 = strfind(pwrstr,']');
% %     if(~isempty(i1))
% %         pwrstr_dest = pwrstr(i1:i2);
% %     else
% %         pwrstr_dest = [];
% %     end
% %     
% %     % add a legend to the dest figure
% %     if(is_hax)
% %         legend(hdestaxes(i),{pwrstr_dest,pwrstr_src})
% %     end
% % end

end
