%function hArray = arrow3d()
    [xcb,ycb]=circle_patch(0,0,1,30);
    zcb = 0*xcb;
    
    hcb = patch(xcb,ycb,zcb,'r');
    
    [xct,yct]=circle_patch(0,0,1,30);
    zct = 0*xct+1;
    
    hct = patch(xct,yct,zct,'r');
    
    patch([xcb(1) xcb xcb(end)],[ycb(1) ycb ycb(end)],[0 0*ycb+1 0],'b');
    
    axis equal;
%end