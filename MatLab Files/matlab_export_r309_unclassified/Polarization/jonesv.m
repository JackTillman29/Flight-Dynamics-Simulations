classdef jonesv < handle
    properties
        A; % unitless fraction of H field magnitude w.r.t. total field
        B; % unitless fraction of V field mantitude w.r.t. total field
        delta; % phase angle of V field w.r.t. H field.
        alpha;
        Eeff;
        Ex;
        Ey;
        deltax;
        deltay;
        J;
        ux,uy,uz;
        az,el;
        S0,S1,S2,S3;
        phi, theta, r;
    end
    methods
        function obj = jonesv(Ex,Ey)
           obj.Ex = Ex;
           obj.Ey = Ey;
           obj.Eeff = sqrt(abs(Ex).^2 + abs(Ey).^2);
           obj.Ex = Ex./obj.Eeff;
           obj.Ey = Ey./obj.Eeff;
           obj.A = abs(Ex)./obj.Eeff;
           obj.B = abs(Ey)./obj.Eeff;
           obj.deltay = angle(Ey);
           obj.deltax = angle(Ex);
           obj.delta  = obj.deltay - obj.deltax;
           obj.alpha = 0.5*atan2((2*obj.A.*obj.B.*cos(obj.delta)),(obj.A.^2-obj.B.^2));
           obj.J = [obj.A obj.B.*exp(1i*obj.delta)].';
           obj.ux = [1 0 0]';
           obj.uy = [0 1 0]';
           obj.uz = [0 0 1]';
           
           % Stokes vector
           %    (x,y,z) on Poincare sphere is (S1,S2,S3)
           %    (S1,S2,S3) is a unit vector
           %     S0 is the magnitude of vector (S1,S2,S3) (S0=1)
           obj.S0 = obj.A.^2 + obj.B.^2;
           obj.S1 = obj.A.^2 - obj.B.^2;
           obj.S2 = 2.*obj.A.*obj.B.*cos(obj.delta);
           obj.S3 = 2.*obj.A.*obj.B.*sin(obj.delta);
               % delta = -chi
           % Stokes vector in spherical coordinates
           [az,el,r] = cart2sph(obj.S1,obj.S2,obj.S3);
           obj.phi   = az;
           obj.theta = el;
           obj.r     = r;
        end
        
        function obj2 = transform(obj,T)
            %gather needed stuff
            new_J = obj.Eeff * T * obj.J;
            new_Eeff = norm(new_J);
            obj2 = jonesv( ...
                exp(1i*angle(obj.Ex)) * new_J(1)./new_Eeff, ...
                                        new_J(2)./new_Eeff); 
        end
        
        function [iq_Ex,iq_Ey] = Polarize(obj,iq)
            % scale to represent net field
%             iq = iq .* obj.Eeff;
            
            % create channels
            iq_Ex = iq .* obj.A;
            iq_Ey = iq .* obj.B .* exp(1i*obj.delta);
        end
        
        function d = Projection(tx,rx,los)
            [d.attenH,d.attenV] = ProjectedElementAreaCompensation(tx,los);
            uxtx = tx.uy;
            uytx = tx.uz;
            uxrx = rx.uy;
            uyrx = rx.uz;
            d.hh = tx.Ex * dot(uxtx,uxrx);
            d.hv = tx.Ey * dot(uytx,uxrx);
            d.vv = tx.Ey * dot(uytx,uyrx);
            d.vh = tx.Ex * dot(uxtx,uyrx);
            d.sx_on_rx_ant = d.attenH * (d.hh + d.hv);
            d.sy_on_rx_ant = d.attenV * (d.vv + d.vh);
            
            d.sx_thru_rx_pol = conj(rx.Ex) * d.sx_on_rx_ant;
            d.sy_thru_rx_pol = conj(rx.Ey) * d.sy_on_rx_ant;
            
            d.gain = abs(d.sx_thru_rx_pol + d.sy_thru_rx_pol).^2;
        end
        
        function [attenH,attenV] = ProjectedElementAreaCompensation(obj,los)
            maglos = sqrt(sum(los.^2));
            uH = obj.uy;
            uV = obj.uz;
            
            attenH = sqrt(sum(cross(uH,los/maglos).^2));
            attenV = sqrt(sum(cross(uV,los/maglos).^2));
        end
        
        function varargout = plot(obj,varargin)
            
            colorspec = 'k';
            if(nargin >= 2)
                if(ishandle(varargin{1}))
                    hfig = varargin{1};
                    newfig = 0;
                    start = 2;
                else
                    start = 1;
                end
                for i = start:2:length(varargin)
                    if(strcmp(varargin{i},'Color'))
                        colorspec = varargin{i+1};
                    else

                    end
                end
                
            else
                newfig = 1;
            end
            
            if(numel(obj.A) == 1)
            
                nPerWave = 100;
                rPerN    = 2*pi./nPerWave;
                nWaves = 5;
                kvec = 0:(nPerWave*nWaves);

                %A represents x,z(k) plane
                x = obj.A * cos(rPerN * kvec);

                %B represents y,z(k) plane rel x
                y = obj.B * cos(rPerN * kvec + obj.delta);

                % quiver x to show DIRECTION of rotation
                karrow = nPerWave/nWaves
                karrowp1 = round(0.1*nPerWave);
                karrowp1 = 1;
                qt = kvec(karrow);
                qx = x(karrow);
                qy = y(karrow);
%                 dt = kvec(karrow+karrowp1) - kvec(karrow);
%                 dx = x(karrow+karrowp1) - x(karrow);
%                 dy = y(karrow+karrowp1) - y(karrow);
                
                dt = kvec(karrow) - kvec(karrow+karrowp1);
                dx = x(karrow) - x(karrow+karrowp1);
                dy = y(karrow) - y(karrow+karrowp1);

    %             obj.alpha = 0.5*atan((2*obj.A*obj.B.*cos(obj.delta))/(obj.A^2-obj.B^2));
    %             obj.alpha = 0.5*atan2((2*obj.A*obj.B.*cos(obj.delta)),(obj.A^2-obj.B^2));
    %             if(isempty(obj.alpha))
    %                 obj.alpha = 0;
    %             end

                if(newfig)
                    figure('Position',[520 148 1000 650]);
                else
                    figure(hfig)
                end
                hax(1) = subplot(1,2,1);
                plot3(kvec,x,0*x,'b:');
                hold on;
                plot3(kvec,0*y,y,'r:');
                plot3(kvec,x,y,'Color',colorspec,'LineWidth',2);
    %             quiver3(qt,qx,qy,dt,dx,dy,'LineWidth',3)
                ht=title(['$$E_{eff}=' num2str(obj.Eeff,'%5.2f') ',A=' num2str(obj.A,'%4.2f') ', B=' num2str(obj.B,'%4.2f') 'e^{j(' num2str(obj.delta) ')}$$'],'Interpreter','latex');
                set(gca,'YDir','reverse');
                hold off;
                xL=xlabel('Propagation Direction \rightarrow');
                set(xL,'Rotation',25);

                ylabel('E_{ox}');
                zlabel('E_{oy}');
                set(gca,'DataAspectRatio',[100 1 1])
                view(-50,20)

                hax(2) = subplot(1,2,2);
                plot3(kvec,x,0*x,'b:');
                hold on;
                plot3(kvec,0*y,y,'r:');
                plot3(kvec,x,y,'Color',colorspec,'LineWidth',2);
                % draw arrow for rotation sense
                ang = 30;
                RT1 = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
                RT2 = [cosd(-ang) -sind(-ang); sind(-ang) cosd(-ang)];
                rotpt1 = RT1*[-dx -dy].';
                rotpt2 = RT2*[-dx -dy].';
                p1 = qx + 2*[0 rotpt1(1)];
                p2 = qy + 2*[0 rotpt1(2)];
                p3 = qx + 2*[0 rotpt2(1)];
                p4 = qy + 2*[0 rotpt2(2)];
                hp1 = plot3([0 0],p1,p2,'k','linewidth',2);
                hp2 = plot3([0 0],p3,p4,'k','linewidth',2);
                ht=title(['$$\alpha=' num2str(180/pi*obj.alpha) '\deg $$'],'Interpreter','latex');
                set(gca,'XDir','normal','YDir','reverse');
%                 view(-90,0);
                view(90,0);
                set(gca,'DataAspectRatio',[500 1 1]);
                ylim([-1 1]);
                zlim([-1 1]);
                ylabel('E_{ox}');
                zlabel('E_{oy}');
                hold off
                try
                    add_print_callbacks;
                catch
                    warning('add_print_callbacks is not in your search path');
                end
            else
                
                
                
            end
            if(nargout > 0)
                varargout{1} = hax;
            end
        end
    end
end