function viz_polar_fft(F,r)
    ang = linspace(-pi,pi,length(F));
    line_base_x = r * cos(ang);
    line_base_y = r * sin(ang);
    line_base_z = 0 * ang;
    
    line_top_z = abs(F);
    
    figure;
    hL=line( ...
        [line_base_x' line_base_x']', ...
        [line_base_y' line_base_y']', ...
        [line_base_z' line_top_z']');
    set(gca,'DataAspectRatio',[1 1 max(abs(F))]);
    set(hL,'Color','b');
end