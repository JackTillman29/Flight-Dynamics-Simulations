function triggerdata = TriggerPoints(filterObj,threshold_on,pretrigger_samples,posttrigger_samples,stream,show)
    if(~exist('show'))
        show = 0;
    end
    
    % for backwards compatibility, theshold_on, if scalar operates the same
    if(length(threshold_on) == 1)
       threshold_on = [threshold_on threshold_on-1];
    end
    
    % arm on the first on->off transition
    if(isa(filterObj,'firFilter'))
        stream_filt = applyFilter(filterObj,stream);
    else
        stream_filt = stream;
    end
    trigBool_rise = stream_filt  > threshold_on(1);
    transitions_rise = max([trigBool_rise 0]-[0 trigBool_rise],0);
    
    trigBool_fall = stream_filt  < threshold_on(2);
    transitions_fall = max([trigBool_fall 0]-[0 trigBool_fall],0);
    
    potentialFall = find(transitions_fall == 1);
    potentialRise = find(transitions_rise == 1);
    
    %transitions_fall(1) = 0;
    %transitions_rise(1) = 0;
    
    % find first rise within each "potential fall" boundary
    potentialTriggers = zeros(1,length(potentialFall)-1);
    for k = 1 : ( length(potentialFall)-1 )
        keep_index = find( ...
        (potentialRise > potentialFall(k)) & ...
        (potentialRise < potentialFall(k+1)) ...
        ,1,'first');
        if(isempty(keep_index))
            potentialTriggers(k) = 0;
        else
            potentialTriggers(k) = potentialRise(keep_index);
        end
    end
    potentialTriggers = potentialTriggers(potentialTriggers ~= 0);
    
    
%     figure;
%     subplot(2,1,1);
%     plot(stream_filt);ah = gca;
%     subplot(2,1,2);
%     plot([trigBool_rise' trigBool_fall']);ah = [ah gca];linkaxes(ah,'x');
%     ylim([-2 2]);
%     return;
    
    

    if(isempty(potentialTriggers))
        triggerdata.lead = [];
        triggerdata.trail = [];
        if(show==1)
            figure;
            plot([stream.' stream_filt.']);
            ylimits = get(gca,'YLim');
            xlimits = get(gca,'XLim');
            line([xlimits],[0 0]+threshold_on(1),'Color','g');
            line([xlimits],[0 0]+threshold_on(2),'Color','r');
            for k = 1 : length(triggerdata.lead)
                line([0 0]+triggerdata.lead(k),ylimits,'Color','g','LineStyle','--');
                line([0 0]+triggerdata.trail(k),ylimits,'Color','r','LineStyle','--');
            end
            return;
        end
    end
    temp = potentialTriggers + posttrigger_samples;
    potentialTriggers = potentialTriggers(temp < length(stream));
    
    temp = potentialTriggers - pretrigger_samples;
    potentialTriggers = potentialTriggers(temp > 0);
    
    
%     if( (potentialTriggers(1) - pretrigger_samples) < 1)
%         potentialTriggers = potentialTriggers(2:end);
%     end
%     if( (potentialTriggers(end)+posttrigger_samples) > length(stream))
%         potentialTriggers = potentialTriggers(1:(end-1));
%     end
    
    triggerdata.lead  = potentialTriggers - pretrigger_samples;
    triggerdata.trail = potentialTriggers + posttrigger_samples;
    
    if(show==1)
    figure;
    %plot([stream.' stream_filt.']);
    plot([stream_filt.']);
    ylimits = get(gca,'YLim');
    xlimits = get(gca,'XLim');
    line([xlimits],[0 0]+threshold_on(1),'Color','g');
    line([xlimits],[0 0]+threshold_on(2),'Color','r');
    for k = 1 : length(triggerdata.lead)
        line([0 0]+triggerdata.lead(k),ylimits,'Color','g','LineStyle','--');
        line([0 0]+triggerdata.trail(k),ylimits,'Color','r','LineStyle','--');
    end
    
end