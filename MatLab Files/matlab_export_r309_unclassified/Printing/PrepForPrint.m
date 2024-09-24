function varargout = PrepForPrint(varargin)
    % Return a default set of properties if nothing passed in
    if(nargin == 0)
        props.titleFontSize   = 16;
        props.xylabelFontSize = 14;
        props.axisFontSize    = 10;
        props.legendFontSize  = 16;
        props.titleFont       = 'Arial';
        props.xylabelFont     = 'Arial';
        props.axisFont        = 'Arial';
        props.legendFont      = 'Arial';
        props.titleFontWgt    = 'bold';
        props.xylabelFontWgt     = 'bold';
        props.axisFontWgt        = 'bold';
        props.legendFontWgt      = 'normal';
        props.lineWidth          = 2;
        props.axesBox            = [.065 .078 1-0.1 1-.129];
        
        props.bgfig = [];
        props.bgax  = [];
        
    % Else use argument 1 as the figure number and argument 2 as the
    % properties structure
        varargout{1} = props;
    else
        % start scanning
        hf = varargin{1};
        props = varargin{2};
        
        %KDS Why filter on empty tags??
        %subPlots = findobj('Type','axes','Tag','','Parent',hf);
        subPlots = findobj('Type','axes','Parent',hf);
        subPlots = [subPlots; findobj('Type','polaraxes','Parent',hf)];
        subPlots = [subPlots; findobj('Type','polaraxes','Parent',hf)];
        disp(['Found ' num2str(length(subPlots)) ' subplots']);
        legends  = findobj('Tag','legend','Parent',hf);
        disp(['Found ' num2str(length(legends)) ' legends']);
        colorbars = findobj(get(hf, 'Children'), 'Tag', 'Colorbar');
        disp(['Found ' num2str(length(colorbars)) ' colorbars']);
        
        % mmpolar is wierd
        mmPolarPlots = findobj('Tag','MMPOLAR_Axes','Parent',hf);
        disp(['Found ' num2str(length(mmPolarPlots)) ' mmpolar axes']);
        
        for k = 1 : length(subPlots)
            ha = subPlots(k);
            hLines = findobj('Type','line','Parent',ha);
            if(~isempty(hLines))
                set(hLines,'LineWidth',props.lineWidth);
            end
            set(ha,'FontSize',props.axisFontSize, ...
                'FontName',props.axisFont, ...
                'FontWeight',props.axisFontWgt);
            try
                hLabels = [get(ha,'XLabel') get(ha,'YLabel') get(ha,'ZLabel')];
            catch
                hLabels = [];
            end
            if(~isempty(hLabels))
                set(hLabels, ...
                'FontSize',props.xylabelFontSize, ...
                'FontName',props.xylabelFont, ...
                'FontWeight',props.xylabelFontWgt);
            end
            set(ha,'FontSize',props.axisFontSize, ...
                'FontName',props.axisFont, ...
                'FontWeight',props.axisFontWgt);
            set( get(ha,'Title'), ...
                'FontSize',props.titleFontSize, ...
                'FontName',props.titleFont, ...
                'FontWeight',props.titleFontWgt);
%             if ( length(subPlots) == 1 )
%                 set(ha,'Position',props.axesBox);
%             end
            
        end
        
        % colorbars
        if(~isempty(colorbars))
            for k = 1:length(colorbars)
                ha = colorbars(k);
                hLabels = [get(ha,'XLabel') get(ha,'YLabel') get(ha,'ZLabel')];
                if(~isempty(hLabels))
                    set(hLabels, ...
                    'FontSize',props.xylabelFontSize, ...
                    'FontName',props.xylabelFont, ...
                    'FontWeight',props.xylabelFontWgt);                    
                end 
            end        
        end
        
        if(~isempty(mmPolarPlots))
            ha = mmPolarPlots;
            hLines = findobj('Type','line','Parent',ha);
            set(hLines,'LineWidth',props.lineWidth);
            hLabels = [get(ha,'XLabel') get(ha,'YLabel') get(ha,'ZLabel')];
%             set(hLabels, ...
%                 'FontSize',props.xylabelFontSize, ...
%                 'FontName',props.xylabelFont, ...
%                 'FontWeight',props.xylabelFontWgt);
%             set(ha,'FontSize',props.axisFontSize, ...
%                 'FontName',props.axisFont, ...
%                 'FontWeight',props.axisFontWgt);
            set( get(ha,'Title'), ...
                'FontSize',props.titleFontSize, ...
                'FontName',props.titleFont, ...
                'FontWeight',props.titleFontWgt);
        end
        
        if(~isempty(legends))
        set(legends,'FontSize',props.legendFontSize, ...
            'FontName',props.legendFont, ...
                'FontWeight',props.legendFontWgt);
        end
        set(hf,'Color','w');
    end
end