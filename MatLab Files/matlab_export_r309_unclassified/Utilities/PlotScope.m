function [hd,hd_ax] = PlotScope(t,dat,...
    currentTime,timeWindowWidth,currentTimePos,...
    datBounds,datStep,...
    datLabel,varargin)

if(length(varargin) > 0)
    iprop = 1;
    for i = 1:2:length(varargin)
        propStr{iprop} = varargin{i};
        propVal{iprop} = varargin{i+1};
        iprop = iprop + 1;
    end
    
end

hd = plot(t, dat, 'b.-');
if(length(varargin) > 0)
    for i = 1:length(propStr)
        set(hd,propStr{i},propVal{i})
    end
end
hd_ax = get(hd,'parent');
try
update_xaxis_with_timestamp(hd_ax)
end
xlim(hd_ax,[currentTime-timeWindowWidth*currentTimePos...
      currentTime+timeWindowWidth*(1-currentTimePos)])
  
if(datBounds(1) ~= 0 & datBounds(2) ~= 0)
    ylim(hd_ax,datBounds)
    set(hd_ax,'ytick',[datBounds(1):datStep:datBounds(2)])
end
grid on
ylabel(datLabel)
                
                
end