% HOW TO USE:
%   script:
%       figure
%       add_group_tool
%       plot(...)
%       ...
%
% Once the figure is ready, you can right click on the edge of the figure
% outside of the axes to access the context menu.
%  1) Select "GroupElements"
%  2) Click and drag on the axes to select the desired elements.
%  3) If the elements are not exactly selected, right click on the blue patch
%     object to modify the width and height and offset the patch
%  4) Once the desired elements are surrounded, right click in the figure
%  outside the axes to access the context menu and select "WriteFile" to
%  get a file that contains:
%       ElementPosX  ElementPosY  GroupNumber
%     A dialog box opens after selecting "WriteFile" for the user to input
%     the filename and GroupNumber.
%  5) If additional groups are desired, then repeat steps 1) through 4) and
%     when you go to write the file, type "0" into the "Create new?" box so
%     that the tool will append the new group to the file instead of
%     overwriting it.


function add_group_tool

hfig = gcf;
% Define a context menu; it is not attached to anything
hcmenu = uicontextmenu;
% Define the context menu items and install their callbacks
item1 = uimenu(hcmenu, 'Label', 'GroupElements', 'Callback', @InteractiveGroupElements);
item2 = uimenu(hcmenu, 'Label', 'WriteFile',     'Callback', @WriteFile)
set(gcf,'uicontextmenu',hcmenu)

set(gcf,'WindowButtonDownFcn',@FigButtonDown)
set(gcf,'WindowButtonUpFcn',@FigButtonUp)
set(gcf,'WindowButtonMotionFcn',@FigButtonMotion)
end







function InteractiveGroupElements(varargin)

set(gca,'NextPlot','add')

xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
w = abs(diff(xlims));
h = abs(diff(ylims));
xini = xlims(1)+w/4;
yini = ylims(1)+h/4;

rect = [xini yini w/2 h/2];
rect_init = [0.25 0.25 0.5 0.5];
set(gcf,'Units','normalized');
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');
temp = rbbox;
point2 = get(gca,'CurrentPoint');

colors = rand(1,3);

% delete old patch object, if there is one
hchild = get(gca,'Children');
if(length(hchild) > 1)
    for ichild = 1:length(hchild)
        if(strcmp(get(hchild(ichild),'Type'),'line'))
            % keep line data
        elseif(~strcmp(get(hchild(ichild),'Tag'),'keep'))
            delete(hchild(ichild))
        end
    end
    hdata = hchild(end);
else
    hdata = hchild(1);
end

rectx = [point1(1,1) point2(1,1) point2(1,1) point1(1,1)];
recty = [point1(1,2) point1(1,2) point2(1,2) point2(1,2)];
midx = (point1(1,1)+point2(1,1))/2;
midy = (point1(1,2)+point2(1,2))/2;
hpatch = patch(rectx,recty,colors,'FaceAlpha',0.5,'ButtonDownFcn',@PatchRightClickCallback)
set(hpatch,'Tag','keep')
groupNumber = inputdlg('Group Number');
groupNumber = str2num(groupNumber{1});
htext = text(midx,midy,num2str(groupNumber));
set(htext,'Tag','keep')
set(hpatch,'UserData',groupNumber);

% find points inside the patch object
minx = min(rectx);
maxx = max(rectx);
miny = min(recty);
maxy = max(recty);

% get data from axes
xdata = get(hdata,'XData');
ydata = get(hdata,'YData');

idx = find( (xdata <= maxx & xdata >= minx) & ...
            (ydata <= maxy & ydata >= miny) );
xsel = xdata(idx);
ysel = ydata(idx);


hp = plot(xsel,ysel,'.','Color',colors);
set(hp,'Tag','keep')

userData = get(gcf,'UserData');
userData.x = xsel;
userData.y = ysel;
userData.group_number(groupNumber) = groupNumber;
userData.group_x{groupNumber} = xsel;
userData.group_y{groupNumber} = ysel;

set(gcf,'UserData',userData)


% % % interactive = 0;
% % % if(interactive)
% % % rect.x = [-1 1 1 -1];
% % % rect.y = [-1 -1 1 1];
% % % rect.top    = max(rect.y);
% % % rect.bottom = min(rect.y);
% % % rect.right  = max(rect.x);
% % % rect.left   = min(rect.x);
% % % rect.centerx = (rect.right + rect.left)/2;
% % % rect.centery = (rect.top   + rect.bottom)/2;
% % % rect
% % % patch(rect.x,rect.y,'b','FaceAlpha',0.5)
% % % 
% % % % draw boxes on top of patch for "drag" functions
% % % scale = 1/100;
% % % prototype.x = [-1 1 1 -1];
% % % prototype.y = [-1 -1 1 1];
% % % % draw top handle
% % % patch(scale*prototype.x+rect.centerx,scale*prototype.y+rect.top,'w','EdgeColor','k')
% % % 
% % % end

set(gca,'NextPlot','replace')

end







function WriteFile(varargin)

userData = get(gcf,'UserData');

if(isfield(userData,'defaultFilename'))
    defaultFilename = userData.defaultFilename;
else
    defaultFilename = '';
end

prompt = {'Filename:';'Group number:';'Create new? (1:yes,0:no):'};
dlg_title = 'Write elements to file...';
num_lines = 1;
defaults = {defaultFilename,'','0'}
answer = inputdlg(prompt,dlg_title,num_lines,defaults);

% write a new file or append to existing
if(str2num(answer{3}))
    fid = fopen(answer{1},'w');
else
    fid = fopen(answer{1},'a');
    userData.defaultFilename = answer{1};
end
groupnum = str2num(answer{2});



% get element locations
x = userData.x;
y = userData.y;

% write locations to file
for i = 1:length(x)
    fprintf(fid,'%d\t%d\t%d\n',x(i),y(i),groupnum);
end

fclose(fid)
set(gcf,'UserData',userData)

end





function point = FigButtonDown(varargin)

point = get(gca,'CurrentPoint');

end



function point = FigButtonUp(varargin)

point = get(gca,'CurrentPoint');

end



function point = FigButtonMotion(varargin)

point = get(gca,'CurrentPoint');

end





function PatchRightClickCallback(varargin)

hobj = varargin{1};
eventData = varargin{2};

% right click is BUTTON 3
if(eventData.Button == 3)
    
    x=get(hobj,'XData');
    y=get(hobj,'YData');
    
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    patchWidth = maxx - minx;
    patchHeight = maxy - miny;
    patchCenter = [(maxx+minx)/2 (maxy-miny)/2];
    
    prompt = {'Width:';'Height:';'Offset X:';'Offset Y:'};
    dlg_title = 'Set width/height of patch object';
    num_lines = 1;
    default_ans = {num2str(patchWidth);num2str(patchHeight);'0';'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,default_ans);
    newWidth = str2num(answer{1});
    newHeight = str2num(answer{2});
    offsetX = str2num(answer{3});
    offsetY = str2num(answer{4});
    
    x = [patchCenter(1)+offsetX-newWidth/2 ...
         patchCenter(1)+offsetX+newWidth/2 ...
         patchCenter(1)+offsetX+newWidth/2 ...
         patchCenter(1)+offsetX-newWidth/2];
    y = [patchCenter(2)+offsetY-newHeight/2 ...
         patchCenter(2)+offsetY-newHeight/2 ...
         patchCenter(2)+offsetY+newHeight/2 ...
         patchCenter(2)+offsetY+newHeight/2];
    set(hobj,'XData',x,'YData',y);
    
elseif(eventData.Button == 1)
    
    groupnum = get(hobj,'UserData');
    defaultFilename = ['group_',num2str(groupnum)];
    prompt = {'Filename:';'Group number:';'Create new? (1:yes,0:no):'};
    dlg_title = 'Write elements to file...';
    num_lines = 1;
    defaults = {defaultFilename, num2str(groupnum), '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaults);

    % write a new file or append to existing
    if(str2num(answer{3}))
        fid = fopen(answer{1},'w');
    else
        fid = fopen(answer{1},'a');
    end
    groupnum = str2num(answer{2});


    userData = get(gcf,'UserData');
    % get element locations
    x = userData.group_x{groupnum};
    y = userData.group_y{groupnum};

    % write locations to file
    for i = 1:length(x)
        fprintf(fid,'%d\t%d\t%d\n',x(i),y(i),groupnum);
    end

    fclose(fid)
    set(gcf,'UserData',userData)
    
end



end




