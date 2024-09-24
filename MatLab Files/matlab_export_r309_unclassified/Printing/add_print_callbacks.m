function add_print_callbacks(varargin)
if(nargin > 0)
    hfig = varargin{1};
else
    hfig = gcf;
end
set(hfig,'Color','w');
% if context menu already exists in figure, get the handle
hcmenu = get(hfig,'uicontextmenu');
% figure doesn't contain a context menu object, create a new one.
%if(~ishandle(hcmenu))
if(isempty(hcmenu))
    hcmenu = uicontextmenu; % not attached to anything yet.
end
% Define callbacks for context menu items that change linestyle
hcb1 = ['if(~exist(''pp'')) pp=PrepForPrint(); end; PrepForPrint(gcf,pp);'];
% Define the context menu items and install their callbacks
item1 = uimenu(hcmenu, 'Label', 'PrepForPrint', 'Callback', hcb1);
item2 = uimenu(hcmenu, 'Label', 'Kprint', 'Callback', 'callback_kprint');
item3 = uimenu(hcmenu, 'Label', 'Snip Tool', 'Callback', 'system(''SnippingTool.exe'');');
item4 = uimenu(hcmenu, 'Label', 'I Feel Lucky', 'Callback', [hcb1 ';callback_kprint(''I Feel Lucky'')']);
set(hfig,'uicontextmenu',hcmenu);

end