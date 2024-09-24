classdef ant2_mod
    
    properties
        
        refFile;
        nAz;
        nEl;
        azPts;
        elPts;
        azGrid; % interp2 format
        elGrid; % interp2 format
        gainDat;
        tilt;
        
    end
    
    methods
        
        function this = ant2_mod(infile)
            
            x = load(infile);
            this.refFile = infile;
            this.nAz = x(1);
            this.nEl = x(2);
            this.elPts = x(3:this.nEl+2);
            rest = x(this.nEl+3:end);
            G = reshape(rest, this.nEl+1, this.nAz);
            this.azPts = G(1,:);
            this.gainDat = G(2:end,:);
            
            %--------------------------------
            if(this.elPts(1) == 0)
                this.elPts = [this.elPts(end:-1:2);this.elPts];
                this.gainDat = [this.gainDat(end:-1:2,:); this.gainDat];
                this.nEl = 2*this.nEl - 1;
                warning('Mirroring Elevation Data')
            end
            %--------------------------------
            
            [this.azGrid,this.elGrid] = meshgrid(this.azPts,this.elPts);
            this.tilt = 0;
            
        end
        
        function plot(this,elcut,azcut)
            
            %elcut is the elevation angle to cut for the azimuth plot
            %azcut is the azimuth angle to cut for the elevation plot
            if(~exist('azcut'))
                azcut = 0;
            end
            if(~exist('elcut'))
                elcut = 0;
            end
            
            figure;
            subplot(1,3,1);
%             imagesc(this.azPts, this.elPts, 10*log10(abs(this.gainDat)));
             surf(this.azPts, this.elPts, 10*log10(abs(this.gainDat)),'EdgeColor','none');
             %imagesc(this.azPts, this.elPts, 10*log10(abs(this.gainDat)));
             view(0,90);
             xlim([-180 180]);ylim([-90 90])
             

            set(gca, 'YDir', 'normal')
            xlabel('Azimuth (deg)')
            ylabel('Elevation (deg)')
            ht = title([this.refFile ' (dBi)']);
            set(ht, 'Interpreter', 'none');
            
            idx = find(this.elPts >= elcut, 1, 'first');
            subplot(1,3,2);
            plot(this.azPts, 10*log10(abs(this.gainDat(idx,:))));
            xlabel('Azimuth (deg)');
            title(['Gain (dBi) @ ' num2str(azcut) '\circ El']);
            grid on;
            
            idx = find(this.azPts >= azcut, 1, 'first');
            subplot(1,3,3);
            plot(this.elPts, 10*log10(abs(this.gainDat(:,idx))));
            xlabel('Elevation (deg)');
            title(['Gain (dBi) @ ' num2str(azcut) '\circ Az']);
            grid on;            
            
            try
                add_print_callbacks
            catch
                disp('Could not add print callbacks! Is it in your path??');
            end
            
        end
        
        function y = getGain(this, AZlookup, ELlookup)
            % returns a gain table equal in size to the input
            % input format should be compatible with
            % [azlookup,ellookup] = meshgrid(...,...)
            
            y = interp2(this.azGrid, this.elGrid+this.tilt, this.gainDat, AZlookup, ELlookup);
            
        end
        
    end
        
end