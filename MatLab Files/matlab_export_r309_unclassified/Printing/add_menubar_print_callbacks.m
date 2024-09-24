function add_menubar_print_callbacks(varargin)
if(nargin > 0)
    hfig = varargin{1};
else
    hfig = gcf;
end
set(hfig,'Color','w');

hmenu = findobj(hfig,'Type','uimenu');
% TODO: probably should do more than check if hmenu is empty because we could
% create other dropdown menus; should eventually check for any of the items
% in this menu defined below.
if(isempty(hmenu))
    % create a menu object
    hmenu = uimenu(hfig,'Label','Printing'); % not attached to anything yet.
    % Define callbacks for context menu items that change linestyle
    hcb1 = ['if(~exist(''pp'')) pp=PrepForPrint(); end; PrepForPrint(gcf,pp);'];
    % Define the context menu items and install their callbacks
    item1 = uimenu(hmenu, 'Label', 'PrepForPrint', 'Callback', hcb1);
    item2 = uimenu(hmenu, 'Label', 'Kprint', 'Callback', 'callback_kprint');
    item3 = uimenu(hmenu, 'Label', 'Snip Tool', 'Callback', 'system(''SnippingTool.exe'');');
    item4 = uimenu(hmenu, 'Label', 'I Feel Lucky', 'Callback', [hcb1 ';callback_kprint(''I Feel Lucky'')']);
end
end