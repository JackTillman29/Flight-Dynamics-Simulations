function append_colorbar_label(hc,str)
    hclab = get(hc,'YTickLabel');
    
    for k = 1 : size(hclab,1)
        hclab2{k} = [hclab(k,:) ' ' str];
    end
    
    set(hc,'YTickLabel',hclab2);
end