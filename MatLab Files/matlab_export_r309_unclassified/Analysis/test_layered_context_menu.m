close all; clear all; clc

figure
plot(randn(1,50),randn(1,50),'.')

% create parent menu
cmenu = uicontextmenu;

% create first category "actions" 
category_actions = uimenu(cmenu,'Label','actions');
action_1 = uimenu(category_actions,'Label','action 1');
action_2 = uimenu(category_actions,'Label','action 2');
action_3 = uimenu(category_actions,'Label','action 3');

% create second category "process"
category_process = uimenu(cmenu,'Label','processes');
process_1 = uimenu(category_process,'Label','process 1');
process_2 = uimenu(category_process,'Label','process 2');

set(gcf,'UIContextMenu',cmenu)





