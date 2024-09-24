function [varargout]= stripchart( ...
    x,xlab,cell_Y,cell_ylab,channelType,titlestr,varargin)
% inputs
%       x: x values for plotting all data contained in cell_Y
% author: Jeff Hole (Booz Allen Hamilton), 2019

if(length(cell_Y) ~= length(cell_ylab) || ...
    length(cell_Y) ~= length(channelType))
    error('cell_Y,cell_ylab, and channelType need to be the same length.')
end

if(length(varargin) ==1)
    subidx = varargin{1};
else
    subidx = 1:length(cell_Y{1});
end

if(length(varargin) ==1)
    yaxisLim = varargin{1};
end

numChannels = length(cell_Y);

%figure;
cnt = 0;
ytick_cnt = 0;
for ichannel = 1:numChannels
    cnt = cnt +1;
    yctr = cnt;
    dat = cell_Y{ichannel};
    if(channelType(ichannel) == 1)
        %Binary Data
        yscale = 0.5;
        bindat = cell_Y{ichannel};
        dat = cell_Y{ichannel}-0.5;
        dat = yscale*dat + yctr;
        %plots dat in gray, and overlays another plot if Only "on" data in 
        %blue
        plot(x,dat,'color',ones(1,3)*0.8)
        hold on; grid off
        tmpx = x;
        tmpy = dat;
        tmpy(bindat == 0) = nan;
        plot(tmpx,tmpy,'b')
        
        %apply tick and tick label information to internal arrays
        ytick_cnt = ytick_cnt + 1;
        ytick(ytick_cnt) = cnt;
        yticklabel{ytick_cnt} = cell_ylab{ichannel};
    
    elseif(channelType(ichannel) == 2)
        %MULTIVALUED ("ANALOG") DATA
        mindat = min(dat(~isinf(dat)));
        maxdat = max(dat(~isinf(dat)));
        if (maxdat > 0 && mindat ~= 0.0)
            maxdat = max(abs(dat));
            mindat = -max(dat);
        end
        extent = maxdat - mindat;
        midpt = (maxdat - mindat)/2;
        if (length(varargin) == 1)
            if (yaxisLim{ichannel} ~= 0.0)
                midpt = 0.0;
                mindat = midpt - yaxisLim{ichannel};
                maxdat = midpt + yaxisLim{ichannel};
                extent = 2 * yaxisLim{ichannel};
            end
        end
        %offset dat so that its min value is 0, then scale so that it's
        % extent is 0.9, then offset to position on strip chart ("yctr")
        
        analog_scale = 0.7;
        dat = analog_scale*(dat - mindat)/extent + yctr - analog_scale/2;
        
        plot(x,dat,'b')
        hold on; grid on
        
        %apply tick and tick label information to internal arrays (INCLUDE
        %LABELS FOR MIN AND MAX OF THE "ANALOG" DATA)
        
        ytick_cnt = ytick_cnt + 1;
        ytick(ytick_cnt) = cnt - analog_scale/2;
        yticklabel{ytick_cnt} = num2str(mindat);
        
        ytick_cnt = ytick_cnt + 1;
        ytick(ytick_cnt) = cnt; 
        yticklabel{ytick_cnt} = cell_ylab{ichannel};
        
        ytick_cnt = ytick_cnt +1;
        ytick(ytick_cnt) = cnt + analog_scale/2;
        yticklabel{ytick_cnt} = num2str(maxdat);
        
    end
end

xlabel(xlab)
title(titlestr, 'Interpreter','none')
%set(gca,'YGrid','on')
set(gca, 'YTick',ytick,'yticklabel',yticklabel)
ylim([0.5 max(ytick)+0.5])

if(nargout == 1)
    varargout{1} = gca;
end

end
