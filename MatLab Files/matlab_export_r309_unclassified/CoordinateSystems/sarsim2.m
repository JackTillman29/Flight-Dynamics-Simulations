% define scatterer locations
span = 1000;
step = 100;
[scatter.x,scatter.y,scatter.z] = meshgrid(-span:step:span,-span:step:span,0);


plat.xi = -45000;
plat.yi = 0*-45000;
plat.zi = 2*9144;
plat.va = 90*pi/180;
plat.vmag  = 1*274;
plat.vz = 0;
plat.vx = plat.vmag * cos(plat.va);
plat.vy = plat.vmag * sin(plat.va);
plat.vel  = [plat.vx plat.vy plat.vz]; 
c = 3.0e8;

radar.pri = 1e-3;
radar.frq = 10e9;
radar.wav = c / radar.frq;
close all;



for t = 0:10:10
    plat.pos = [ ...
        plat.xi + t*plat.vx ...
        plat.yi + t*plat.vy ...
        plat.zi + t*plat.vz ];

    los.x = scatter.x - plat.pos(1);
    los.y = scatter.y - plat.pos(2);
    los.z = scatter.z - plat.pos(3);
 
    los.xc = -plat.pos(1); 
    los.yc = -plat.pos(2);
    los.zc = -plat.pos(3);
    
    M = SAR_PointingVectors( [los.xc los.yc los.zc]./norm(plat.pos));
    if(t == 0)
        figure(99);
        h=M*[los.x(:) los.y(:) los.z(:)]';
        
        % h represents los,crossrng,downish
        h(3,:) = 0;
        h2 = M'*h;
        % h2 is back in inertial with no height
        
        plot3(plat.pos(1)+h2(1,:),plat.pos(2)+h2(2,:),plat.pos(3)+h2(3,:),'.');
        
        hold on;
        plot3(scatter.x(:),scatter.y(:),scatter.z(:),'r.');
        %line([plat.xi 0],[plat.yi 0],[plat.zi 0]);
        
        cf = CartesianFrame('SAR');
        cf.Rotate(M);
        draw(cf,[0 0 0],1000);
        axis equal;
    
        cf_c = CartesianFrame('PC');
        draw(cf_c,[0 0 0],1000);
    end

    % M = [ crossrange | LOS | normal ]
    
    plat.vmag_tan =  dot(plat.vel,M(2,:));

    scatter.range = sqrt( ...
        (los.x).^2 + ...
        (los.y).^2 + ...
        (los.z).^2 );

    scatter.range_c = sqrt( ...
        (los.xc).^2 + ...
        (los.yc).^2 + ...
        (los.zc).^2 );
    scatter.range_ground_c = sqrt( ...
        (los.xc).^2 + ...
        (los.yc).^2 );

    signal.toa =   2*scatter.range   ./ c;
    signal.toa_c = 2*scatter.range_c ./ c;
    
    signal.relvel = (los.x * plat.vx + los.y * plat.vy + los.z * plat.vz) ./ ...
        scatter.range;
    
    signal.relvel_c = (los.xc * plat.vx + los.yc * plat.vy + los.zc * plat.vz) ./ ...
        scatter.range_c;
    
    signal.doppler = 2.0 * (signal.relvel-signal.relvel_c) / radar.wav;
    
    signal.downrange  = (c*signal.toa./2);
    signal.groundrange = sqrt( signal.downrange.^2 - plat.pos(3).^2);  
    signal.crossrange = signal.downrange .* asin(radar.wav.*signal.doppler./(2*plat.vmag));
    
    
    signal.downrange_proj = sqrt( signal.groundrange.^2 - signal.crossrange.^2 ) - scatter.range_ground_c; 

    plat.zangrate = plat.vmag_tan./scatter.range_c;
    
    % rotate x,y based on scene and platform velocity
    rA = plat.zangrate * t;
    %rA = 0;
    ca = cos(rA);
    sa = sin(rA);
    RM = [ca sa;-sa ca];
    
    tempx = signal.crossrange(:);
    tempy = signal.downrange(:);
    temp = RM * [tempx tempy]';
    signal.downrange_rot = reshape(temp(2,:),size(signal.downrange,1),size(signal.downrange,2));
    
    
    figure(1);
    subplot(1,2,1);
    plot(signal.doppler(:),signal.downrange(:),'.');
    hold on;
    
    subplot(1,2,2);
    plot(signal.crossrange(:),signal.downrange(:),'.');
    hold on;
    plot(signal.crossrange_rot(:),signal.downrange_proj_rot(:),'o');
    drawnow;
    
    figure(2);
    plot(t,plat.vmag_tan,'b.');
    hold on;
    
    figure(3);
    plot(t,plat.zangrate,'b.');
    hold on;
    
end
figure(1);
subplot(1,2,1);
xlabel('Doppler (Hz)');
ylabel('Range (m)');
subplot(1,2,2);
xlabel('Crossrange (m)');
ylabel('Downrange (m)');

figure(2);
xlabel('Time (s)');
ylabel('Tangential Velocity (m/s)');

figure(3);
xlabel('Time (s)');
ylabel('Angular Rate (rad/s)');
