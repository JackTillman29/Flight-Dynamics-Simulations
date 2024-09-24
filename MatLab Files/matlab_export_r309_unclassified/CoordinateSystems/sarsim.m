% define scatterer locations
close all;clear all;clc;
span =1000;
step = 100;
[scatter.x,scatter.y,scatter.z] = meshgrid(-span:step:span,-span:step:span,0);

filedat = load('z:\sawmillkd\triangle2.dat');
scatter.x = filedat(:,1);
scatter.y = filedat(:,2);
scatter.z = filedat(:,3);


plat.xi = 45000;
plat.yi = 45000;
plat.zi = 1*9144;
plat.va = (90+45)*pi/180;
plat.vmag  = 1*300;
plat.vz = 0;
plat.vx = plat.vmag * cos(plat.va);
plat.vy = plat.vmag * sin(plat.va);
plat.vel  = [plat.vx plat.vy plat.vz]; 
c = 3.0e8;

radar.pri = 1e-3;
radar.frq = 10e9;
radar.wav = c / radar.frq;
close all;
dt = 1;
time_vec = 0:dt:32;

omega_vec = 0*time_vec;
k = 1;
for t = time_vec
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
        temp = [los.x(:) los.y(:) los.z(:)];
        temp_rng = sqrt(sum(temp.^2,2));
        h=M*([ ...
            los.x(:) ...
            los.y(:) ...
            los.z(:)])';
        
        % because the SAR frame has no ability to measure "down"
        h(3,:) = 0;
        
        % generate unit vectors to each scatterer in SAR frame
        for ii = 1 : size(h,2)
            h(:,ii) = h(:,ii) ./ norm(h(:,ii));
        end
        
        % convert unit vectors to inertial frame
        h2 = (M'*h)';
        
        los.xsar = h2(:,1) .* temp_rng;
        los.ysar = h2(:,2) .* temp_rng;
        los.zsar = h2(:,3) .* temp_rng;
        
        
        % h2 is back in inertial with no height
        plot3( ...
            plat.pos(1)+los.xsar, ...
            plat.pos(2)+los.ysar, ...
            plat.pos(3)+los.zsar, ...
            'r.');
        
        hold on;
        plot3(plat.pos(1)+los.x,plat.pos(2)+los.y,plat.pos(3)+los.z,'b.');
        %line([plat.xi 0],[plat.yi 0],[plat.zi 0]);
        
        cf = CartesianFrame('SAR',M);
        draw(cf,M'*[-200 0 0]',100);
        axis equal;
    
         cf_c = CartesianFrame('PC');
         draw(cf_c,[0 0 0],100);
    end
 
    % M = [ downrange; crossrange; downish ]
    
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
    
    signal.doppler = 2.0 * signal.relvel / radar.wav;
    signal.doppler_c = 2.0 * signal.relvel_c / radar.wav;
    signal.range = 0.5*c*signal.toa;
    
    signal.groundrange = sqrt( signal.range.^2 - plat.pos(3).^2);
   
    signal.crossrange = signal.range .* asin(radar.wav.*(signal.doppler-signal.doppler_c)./(2*plat.vmag));
    signal.downrange_proj  = sqrt( signal.groundrange.^2 - signal.crossrange.^2 ); 
    
    % make relative to scene center
    signal.downrange_proj = signal.downrange_proj - scatter.range_ground_c;

    omega_vec(k) = plat.vmag_tan./scatter.range_ground_c;
   
    % rotate x,y based on scene and platform velocity
    rA = sum(omega_vec(1:(k-1)).*dt);
    %rA = omega_vec(k) * t;
    %rA = -0.06;
    %rA = 0;
    ca = cos(rA);
    sa = sin(rA);
    RM = [ca -sa;sa ca];
    
    tempx = signal.crossrange(:);
    tempy = signal.downrange_proj(:);
    temp = RM * [tempx tempy]';
    signal.downrange_proj_rot = reshape(temp(2,:),size(signal.downrange_proj,1),size(signal.downrange_proj,2));
    signal.crossrange_rot = reshape(temp(1,:),size(signal.crossrange,1),size(signal.crossrange,2));
    
    figure(1);
    subplot(2,2,[1 3]);
    plot(signal.doppler(:)-signal.doppler_c,signal.downrange_proj(:),'.');
    hold on;
    
    subplot(2,2,2);
    plot(signal.crossrange(:),signal.downrange_proj(:),'.');
    hold on;axis equal;title('Uncompensated Downrange vs Crossrange');
    subplot(2,2,4);
    plot(signal.crossrange_rot(:),signal.downrange_proj_rot(:),'rx');
    axis equal;hold on;title('Compensated Downrange vs Crossrange');
    drawnow;
    
    %figure(2);
    %plot(t,plat.vmag_tan,'b.');
    %hold on;
    
    %figure(3);
    %plot(t,omega_vec(k),'b.');
    %hold on;
    k = k + 1;
end
figure(1);
subplot(2,2,[1 3]);
xlabel('Doppler (Hz)');
ylabel('Range (m)');
subplot(2,2,2);
xlabel('Crossrange (m)');
ylabel('Downrange (m)');
axis equal;
subplot(2,2,4);
xlabel('Crossrange (m)');
ylabel('Downrange (m)');
axis equal;

% figure(2);
% xlabel('Time (s)');
% ylabel('Tangential Velocity (m/s)');
% 
% figure(3);
% xlabel('Time (s)');
% ylabel('Angular Rate (rad/s)');
