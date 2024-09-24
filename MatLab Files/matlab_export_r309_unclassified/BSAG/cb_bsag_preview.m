function varargout = cb_bsag_preview(varargin)
 k = waitforbuttonpress;
 point1 = get(gca,'CurrentPoint');    % button down detected
 finalRect = rbbox;                   % return figure units
 point2 = get(gca,'CurrentPoint');    % button up detected
 idx_start = point1(1,1);
 idx_span  = point2(2,1) - idx_start;
end