function cm = colormap_zeroblack(N,posColor,negColor)


poscm = [...
    linspace(posColor(1),0,floor(N/2)).' ...
    linspace(posColor(2),0,floor(N/2)).' ...
    linspace(posColor(3),0,floor(N/2)).'];

negcm = [...
    linspace(0,negColor(1),floor(N/2)).' ...
    linspace(0,negColor(2),floor(N/2)).' ...
    linspace(0,negColor(3),floor(N/2)).'];

cm = [poscm; negcm];
cm = cm(end:-1:1,:);


end