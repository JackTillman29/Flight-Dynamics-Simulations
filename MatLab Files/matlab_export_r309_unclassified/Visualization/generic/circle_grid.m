function circle_grid(str)
loc_xtick = get(gca,'XTick');
loc_ytick = get(gca,'YTick');

if(mean(diff(loc_xtick)) > mean(diff(loc_ytick)))
    % use ytick
    ref_tick = loc_ytick;
else
    % use xtick
    ref_tick = loc_xtick;
end

ref_tick = ref_tick(ref_tick > 0);
dref = mean(diff(ref_tick));
ref_tick = [ ref_tick ref_tick(end)+dref ref_tick(end)+2*dref ];

hpc = [];
hpcsub = [];
ht_label = [];
for k = 1 : length(ref_tick)
    [xc,yc] = circle_patch(0,0,ref_tick(k),50);
    hpc = [hpc patch(xc,yc,'r')];
    if(k ~= length(ref_tick))
    ht_label = [ ht_label text( ...
        ref_tick(k)*cos(45*pi/180), ...
        ref_tick(k)*sin(45*pi/180), ...
        [num2str(ref_tick(k)) str])];
    % half circle
    span = ref_tick(k+1)-ref_tick(k);
    [xc,yc] = circle_patch(0,0,ref_tick(k)+span/2,50);
    hpcsub = [hpcsub patch(xc,yc,'r')];
    end
end
set(hpc,'FaceColor','none','LineStyle','-','EdgeColor',[0 0 0]+0.5);
set(hpcsub,'FaceColor','none','LineStyle',':','EdgeColor',[0 0 0]+0.5);
set(ht_label,'Color',[0 0 0]+1.0, 'FontSize', 12, 'FontWeight', 'bold');
set([hpc hpcsub ht_label],'Tag','circle_grid');
end