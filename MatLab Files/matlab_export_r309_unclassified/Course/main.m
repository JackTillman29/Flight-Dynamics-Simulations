close all; clear all; clc


% vidobj = VideoWriter('video.mp4');

make_gif = 1;
    gif_filename = 'test.gif';


% unit circle
N = 2e3;
th = linspace(-pi,pi,N);
ci = cos(th);
cq = sin(th);

% sinusoid
Fs = 100; dt = 1/Fs;
f = 2;
nframes = 100;
t = [0:(nframes-1)]*dt;

posf = exp(1i*2*pi*f.*t);
negf = exp(-1i*2*pi*f.*t);
sumf = 0.5*(posf + negf);

colors = lines;
gray_color = 0.8*ones(1,3);


% fmt_str = '%10.1f';
% str_function = @(x) 

for iframe = 1:nframes

    if(iframe == 1)
        hfig = figure('Position',[237 410 1485 420]);
        set(gcf,'color','w')
        subplot(1,3,1);
        plot(ci,cq,'color',gray_color);
        set(gca,'xcolor','w','ycolor','w')
        hold on; axis equal;
        plot([0 0],[0 1.25],'--','color',gray_color); text(0.01,1.2,'Q','Color',gray_color,'FontWeight','bold','FontSize',12);
        plot([0 1.25],[0 0],'--','color',gray_color); text(1.2,0.1,'I','Color',gray_color,'FontWeight','bold','FontSize',12);
        
        hp(1) = plot(real(posf(iframe)),imag(posf(iframe)),'Marker','o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
        hp(2) = plot([0 real(posf(iframe))],[0 imag(posf(iframe))],'Color',colors(1,:));
        htext(1) = text(0,-0.5,str_function(posf(iframe)),'FontSize',14,'FontWeight','bold','FontAngle','italic','HorizontalAlignment','center','Color',gray_color);
        xlim(1.25*[-1 1])
        ylim(1.25*[-1 1])
        
        subplot(1,3,2);
        plot(ci,cq,'color',gray_color);
        set(gca,'xcolor','w','ycolor','w')
        hold on;
        plot([0 0],[0 1.25],'--','color',gray_color); text(0.01,1.2,'Q','Color',gray_color,'FontWeight','bold','FontSize',12);
        plot([0 1.25],[0 0],'--','color',gray_color); text(1.2,0.1,'I','Color',gray_color,'FontWeight','bold','FontSize',12);
        
        hp(3) = plot(real(negf(iframe)),imag(negf(iframe)),'Marker','o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
        hp(4) = plot([0 real(negf(iframe))],[0 imag(negf(iframe))],'Color',colors(2,:));
        htext(2) = text(0,-0.5,str_function(negf(iframe)),'FontSize',14,'FontWeight','bold','FontAngle','italic','HorizontalAlignment','center','Color',gray_color);
        
        xlim(1.25*[-1 1])
        ylim(1.25*[-1 1])
        
        text(-1.65,0,'+','FontSize',64,'FontWeight','bold','Color',gray_color,'HorizontalAlignment','center')
        
        
        subplot(1,3,3);
        plot(ci,cq,'color',gray_color);
        set(gca,'xcolor','w','ycolor','w')
        hold on;
        plot([0 0],[0 1.25],'--','color',gray_color); text(0.01,1.2,'Q','Color',gray_color,'FontWeight','bold','FontSize',12);
        plot([0 1.25],[0 0],'--','color',gray_color); text(1.2,0.1,'I','Color',gray_color,'FontWeight','bold','FontSize',12);
        
        hp(5) = plot(real(posf(iframe)),imag(posf(iframe)),'Marker','o','MarkerFaceColor',gray_color,'MarkerEdgeColor','none');
        hp(6) = plot(real(negf(iframe)),imag(negf(iframe)),'Marker','o','MarkerFaceColor',gray_color,'MarkerEdgeColor','none');
        hp(7) = plot(real(sumf(iframe)),imag(sumf(iframe)),'Marker','o','MarkerFaceColor',zeros(1,3),'MarkerEdgeColor','none');
        
        htext(3) = text(0,-0.5,str_function(sumf(iframe)),'FontSize',14,'FontWeight','bold','FontAngle','italic','HorizontalAlignment','center','Color',gray_color);
        
        hp(8) = plot([0 real(posf(iframe))],[0 imag(posf(iframe))],'Color',gray_color);
        hp(9) = plot([0 real(negf(iframe))],[0 imag(negf(iframe))],'Color',gray_color);
        hp(10) = plot([0 real(sumf(iframe))],[0 imag(sumf(iframe))],'Color',zeros(1,3));
        xlim(1.25*[-1 1])
        ylim(1.25*[-1 1])
        
        text(-1.65,0,'=','FontSize',64,'FontWeight','bold','Color',gray_color,'HorizontalAlignment','center')
        
        
    else
        set(hp(1),'XData',real(posf(iframe)),'YData',imag(posf(iframe)))
        set(hp(2),'XData',[0 real(posf(iframe))],'YData',[0 imag(posf(iframe))])
        
        set(hp(3),'XData',real(negf(iframe)),'YData',imag(negf(iframe)))
        set(hp(4),'XData',[0 real(negf(iframe))],'YData',[0 imag(negf(iframe))])
        
        set(hp(5),'XData',real(posf(iframe)),'YData',imag(posf(iframe)))
        set(hp(6),'XData',real(negf(iframe)),'YData',imag(negf(iframe)))
        set(hp(7),'XData',real(sumf(iframe)),'YData',imag(sumf(iframe)))
        
        set(hp(8),'XData',[0 real(posf(iframe))],'YData',[0 imag(posf(iframe))])
        set(hp(9),'XData',[0 real(negf(iframe))],'YData',[0 imag(negf(iframe))])
        set(hp(10),'XData',[0 real(sumf(iframe))],'YData',[0 imag(sumf(iframe))])
        
        set(htext(1),'String',str_function(posf(iframe)))
        set(htext(2),'String',str_function(negf(iframe)))
        set(htext(3),'String',str_function(sumf(iframe)))
    end
    
    if(make_gif)
        
        this_frame = getframe(hfig);
        im = frame2im(this_frame);
        [imind,cm] = rgb2ind(im,256);
        
        if( iframe == 1 )
            imwrite(imind,cm,gif_filename,'gif', 'LoopCount', inf);
        else
            imwrite(imind,cm,gif_filename,'gif', 'WriteMode', 'append');
        end
    end
    
    drawnow;
    pause(0.1)
%     return

end





