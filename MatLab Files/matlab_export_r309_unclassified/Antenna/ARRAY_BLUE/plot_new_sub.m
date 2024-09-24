nSubs = 16;


figure;
hold on;
for ksub = 1 : nSubs
    yd = eval(sprintf('new_sub%d(:,2);',ksub));
    zd = eval(sprintf('new_sub%d(:,3);',ksub));
    hp=plot(yd,zd,'.','MarkerSize',18);
    disp(length(zd))
end
set(gca,'DataAspectRatio',[1 1 1]);