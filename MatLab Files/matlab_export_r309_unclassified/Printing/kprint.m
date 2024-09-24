function kprint(fignum,filename,dpi)
if(~isnumeric(fignum))
    fignum = fignum.Number;
end
set(fignum,'PaperPositionMode','auto');
% NEED TO TURN OFF 'InvertHardCopy' PROPERTY IN THE FIGURE PROPERTIES,
% OTHERWISE THE PRINT COMMAND WILL IGNORE AXES BACKGROUND COLORS AND PRINT
% ONLY WHITE BACKGROUNDS
set(fignum,'InvertHardCopy','off');

bo = 5; % pixel border
disp(['Creating PNG baseline file for ' filename]);

% ensure temp file is unique (no pre-existing __temp.png)
%    had issues using kprint with multiple instances of matlab running in
%    parallel using the same working directory
% tmpfilename = '__temp.png';
tmpfilename = ['__temp_',num2str(rand,'%10.8f'),'.png'];
% itmp = 0;
% while(exist(tmpfilename,'file'))
%     itmp = itmp + 1;
%     if(~exist(tmpfilename,'file'))
%         break
%     else
%         % increment tmp number and try this filename
%         tmpfilename = ['__temp',num2str(itmp),'.png'];
%     end
% end

print('-dpng',['-f' num2str(fignum)],['-r' num2str(dpi)],tmpfilename);
disp(['Post-Processing to remove exterior white space']);
q2 = imread(tmpfilename);
% q2 = getframe(fignum);
% q2 = q2.cdata;
qsum = sum(q2,3);
[xi,xj]=find(qsum~=765);
try
size_q2 = size(q2);
qmod = q2( ...
    max((min(xi)-bo),1) : min((max(xi)+bo), size_q2(1)), ...
    max((min(xj)-bo),1) : min((max(xj)+bo), size_q2(2)), :);
q2 = qmod;
catch
    disp('Issues removing whitespace. Proceeding unaltered.');
end

imwrite(q2,filename,'png');
delete(tmpfilename);
end