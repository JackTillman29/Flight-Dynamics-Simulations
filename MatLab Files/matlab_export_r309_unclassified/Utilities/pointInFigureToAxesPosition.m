function pt_in_fig = pointInFigureToAxesPosition(pt_in_fig,hax)
% author: Jeff Hole (Booz Allen Hamilton)
% 11-19-2019

axes_pos_in_fig = get(hax,'position');

axes_xlim = get(hax,'xlim');
axes_ylim = get(hax,'ylim');
axes_pos_in_axes(1) = axes_xlim(1);
axes_pos_in_axes(2) = axes_ylim(1);
axes_pos_in_axes(3) = abs(diff(axes_xlim));
axes_pos_in_axes(4) = abs(diff(axes_ylim));

% scale and shift
pt_in_axes_norm(:,1) = (pt_in_fig(:,1) - axes_pos_in_fig(1)) ./ axes_pos_in_fig(3);
pt_in_axes_norm(:,2) = (pt_in_fig(:,2) - axes_pos_in_fig(2)) ./ axes_pos_in_fig(4);
pt_in_fig(:,1) = pt_in_axes_norm(:,1) .* axes_pos_in_axes(3) + axes_pos_in_axes(1);
pt_in_fig(:,2) = pt_in_axes_norm(:,2) .* axes_pos_in_axes(4) + axes_pos_in_axes(2);






end