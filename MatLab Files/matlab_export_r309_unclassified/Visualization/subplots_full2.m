function handles_return = subplots_full2(figh,nh,nv,b)
    % jah: added ability for specifying figure handle, and needed to flip 
    % bottom nested loop v and h vectors for alignment...
    
    if(isnumeric(figh))
        if(figh ~= 0)
            figure(figh)
        end
    else
        figure(figh)
    end
    
    nsph = nh; % num horizontal
    nspv = nv; % num vertical
    fsph = 1/nsph; % fraction of fig horz
    fspv = 1/nspv; % fraction of fig vert
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
                0+(kv-1)*fspv+0.5*b*fspv ...   % startx ( flip h -> v
                1-kh*fsph+0.5*b*fsph ...       % starty ( flip v -> h
                fspv-b*fspv ...            % width      ( flip h -> v
                fsph-b*fsph]);             % height     ( flip v -> h
            iwin = iwin + 1;
        end
    end
    handles_return = ah;
end