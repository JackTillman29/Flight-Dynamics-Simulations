function linkallaxes(figlist,linkarg)
    ah = [];
    % build a list of all axes handles associated with the figure #s passed in
    for k = 1 : length(figlist)
        ah = [ah findobj('parent',figlist(k),'Type','axes')];
    end
    
    % call linkaxes as usual
    linkaxes(ah,linkarg);
end


