function handles_return = subplots_full(figh,nh,nv,b)
    nsph = nh; % num horizontal
    nspv = nv; % num vertical
    fsph = 1/nsph; % fraction of fig horz
    fspv = 1/nspv; % fraction of fig vert
    % b is border
    iwin = 1;
    ah = zeros(nsph,nspv);
    for kh = 1 : nsph
        for kv = 1 : nspv
            ah(kh,kv)=subplot(nsph,nspv,iwin);
            set(gca,'ActivePositionProperty','OuterPosition');
            iwin = iwin + 1;
        end
    end

    iwin = 1;
    for kh = 1 : nsph
        for kv = 1 : nspv
            set(ah(kh,kv),'Visible','off','Position',[ ...
                0+(kh-1)*fsph+0.5*b*fsph ...   % startx
                1-kv*fspv+0.5*b*fspv ...       % starty
                fsph-b*fsph ...            % width
                fspv-b*fspv]);             % height
            iwin = iwin + 1;
        end
    end
    handles_return = ah;
end