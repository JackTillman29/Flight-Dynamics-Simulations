function get_on_screen_power( varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    uih=get(gcbo,'Parent');
    udat = get(uih,'UserData');
    ah = udat{1};
    Lh = get(ah,'Children');
    totpow = udat{2};
    xlimits = get(ah,'XLim');
    xd = get(Lh,'XData');
    yd = get(Lh,'YData');
    idx1 = find(xd > xlimits(1),1,'first');
    idx2 = find(xd < xlimits(2),1,'last');
    thispow = totpow*(yd(idx2)-yd(idx1))

end

