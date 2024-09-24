    clear;
    clc;
    close all;
    
    c = 3e8;
    r2d = 180/pi;
    d2r = pi/180;
    
    % define antenna
    load ant.mat;
    figure;imagesc(ant.gain);
    caxis([-40 0])
    return
    rf = 35e9;
    % wavelen = 0.03;
    wavelen = c/rf;
    ddop = 500;
    
    % ned
    % mis.pos = [0 0 12192];
    mis.pos = [0 0 9144];
    % mis.vel = [274 0 0];
    % mis.vel = [1500 0 0];
    mis.speed = 1500;
    mis.yaw   = 0*pi/180;
    mis.pitch = -30*pi/180;
    
    % yaw = linspace(0,30,30);
    % pitch = linspace(-30,-10,30);
    % for j = 1:30
    %     
    %     mis.yaw = yaw(j) * d2r;
    %     mis.pitch = pitch(j) * d2r;
    
    mis.vel = mis.speed*[cos(mis.yaw)*cos(mis.pitch) sin(mis.yaw)*cos(mis.pitch) sin(mis.pitch)];
    
    % ASSUMING LOS == MISSILE BODY FRAME (VELOCITY VECTOR, defined by pitch/yaw)
    
    % Define aperture beam
    mis.azbw = 1 * pi/180;
    mis.elbw = 1 * pi/180;
    
    % Define aperture LOS vectors
    nLOS = 8; % # of edge-of-beam LOS vectors
    mis.los_vectors = zeros(1+nLOS,3);
    [x,y,z]= sph2cart(mis.yaw,mis.pitch,1);
    mis.los_vectors(1,:) = [x y z];
    
    
    th = linspace(0,2*pi,nLOS+1); th = th(1:end-1);
    for klos = 1:nLOS
        [x,y,z]= sph2cart( ...
        mis.yaw+sin(th(klos))*mis.azbw/2, ...
        mis.pitch+cos(th(klos))*mis.elbw/2,1);
        mis.los_vectors(1+klos,:) = [x y z];
    end
    
    
    % [x,y,z]= sph2cart( ...
    %     mis.yaw+mis.azbw/2, ...
    %     mis.pitch+mis.elbw/2,1);
    % mis.los_vectors(2,:) = [x y z];
    % 
    % [x,y,z]= sph2cart( ...
    %     mis.yaw+mis.azbw/2, ...
    %     mis.pitch-mis.elbw/2,1);
    % mis.los_vectors(3,:) = [x y z];
    % 
    % [x,y,z]= sph2cart( ...
    %     mis.yaw-mis.azbw/2, ...
    %     mis.pitch+mis.elbw/2,1);
    % mis.los_vectors(4,:) = [x y z];
    % 
    % [x,y,z]= sph2cart( ...
    %     mis.yaw-mis.azbw/2, ...
    %     mis.pitch-mis.elbw/2,1);
    % mis.los_vectors(5,:) = [x y z];
    clear x y z;
    
    gndIntersectPt = zeros(nLOS+1,3);
    for k = 1:nLOS+1
        gndIntersectPt(k,1:3) = mis.pos - mis.pos(3) / mis.los_vectors(k,3) * mis.los_vectors(k,:);
    end
    
    %%
    % figure; plot3(mis.pos(1),mis.pos(2),mis.pos(3),'k.','MarkerSize',15);
    % hold on;
    % plot3([mis.pos(1)],mis.pos(2),mis.pos(3),'k.','MarkerSize',15);
    
    dr = 100;
    % custom grid==============================================================
    % xa = 0:dr:40e3;
    % ya = -20e3:dr:20e3;
    % grid centered on ground intersection point===============================
    gridextent = 20e3;
    xa = gndIntersectPt(1,1) + [-gridextent:dr:gridextent];
    ya = gndIntersectPt(1,2) + [-gridextent:dr:gridextent];
    N = length(xa) * length(ya);
    [xg,yg]=meshgrid(xa,ya);
    
    % compute missile to ground
    vec_m2g_x = xg-mis.pos(1);
    vec_m2g_y = yg-mis.pos(2);
    vec_m2g_z = 0 -mis.pos(3);
    
    
    [azg,elg,rg]=cart2sph(vec_m2g_x,vec_m2g_y,vec_m2g_z);
    % add missile pitch and yaw
    azg = azg-mis.yaw;
    elg = elg-mis.pitch;
    clutter.rg = rg;
    maxr = max(rg(:));
    
    clutter.ratten = -20*log10(4*pi*rg.*rg);
    clutter.antgain = (interp2(ant.el,ant.az,ant.gain,elg,azg))+30;
    clutter.power = clutter.ratten + clutter.antgain;
    clutter.antgain_peak = max(max(clutter.antgain));
    
    if(j==1)
    figure; add_print_callbacks;
    end
    imagesc(1e-3*xg(1,:),1e-3*yg(:,1),clutter.antgain);
    hold on;
    for k = 1:(nLOS+1)
        plot(1e-3*gndIntersectPt(k,1),1e-3*gndIntersectPt(k,2),'k.','MarkerSize',15)
    end
    hold off;
    caxis([-10 30]);
    xlabel('DR (km)');ylabel('CR (km)');
    title('Projected Antenna Gain on Ground (dBi)');
    
    drawnow;
    
    % now compute doppler for each cell
    cvel = -(vec_m2g_x * mis.vel(1) + vec_m2g_y * mis.vel(2) + vec_m2g_z * mis.vel(3)) ./ rg;
    clutter.dop  = -2*cvel ./ wavelen;
    
    figure; add_print_callbacks;
    imagesc(1e-3*xg(1,:),1e-3*yg(:,1),1e-3*clutter.dop);
    hold on
    [junk,h] = contour(1e-3*xg,1e-3*yg,1e-3*clutter.dop,10);
    set(h,'LineColor','k')
    xlabel('DR (km)');ylabel('CR (km)');
    title('Doppler Return from Ground (kHz)');colorbar;
    
    figure; add_print_callbacks;
    imagesc(1e-3*xg(1,:),1e-3*yg(:,1),1e-3*clutter.rg);
    xlabel('DR (km)');ylabel('CR (km)');
    title('Range of Return from Ground (km)');colorbar;
    
    rbin = 0:500:max(clutter.rg(:));
    dbin = 0:500:max(clutter.dop(:));
    bindata = zeros(length(rbin),length(dbin)) * nan;
    
    %% side cut
    figure; add_print_callbacks;
    plot(1e-3*clutter.dop(101,:),clutter.power(101,:));
    grid on;
    title('Clutter Doppler - X Axis Only');
    xlabel('kHz');
    ylabel('dB');
    
    figure; add_print_callbacks;
    plot(1e-3*clutter.rg(101,:),clutter.power(101,:));
    grid on;
    title('Clutter Range - X Axis Only');
    xlabel('km');
    ylabel('dB');
    
    figure; add_print_callbacks;
    plot(1e-3*clutter.rg(101,:),1e-3*clutter.dop(101,:));
    grid on;
    title('Clutter Doppler vs Range - X Axis Only');
    xlabel('km');
    ylabel('kHz');
    
    drawnow;
    
    %%
    
    % % % % % disp('1st method...')
    % % % % % tic
    % % % % % numFound = zeros(length(rbin),length(dbin));
    % % % % % for kr = 1 : length(rbin) - 1
    % % % % %     for kd = 1 : length(dbin) - 1
    % % % % %         res = find( ...
    % % % % %             clutter.rg > rbin(kr) & ...
    % % % % %             clutter.rg < rbin(kr+1) & ...
    % % % % %             clutter.dop > dbin(kd) & ...
    % % % % %             clutter.dop < dbin(kd+1));
    % % % % %         if(~isempty(res))
    % % % % %             numFound(kr,kd) = length(res);
    % % % % %             % get maximum for now
    % % % % %             bindata(kr,kd) = max(clutter.power(res));
    % % % % %         end
    % % % % %     end
    % % % % %     %disp(kr./length(rbin));
    % % % % % end
    % % % % % toc
    % % % % % figure; add_print_callbacks;
    % % % % % imagesc(1e-3*dbin,1e-3*rbin,bindata);
    % % % % % xlabel('Unambiguous Doppler (kHz)');
    % % % % % ylabel('Unambiguous Range (km)');
    % % % % % title('Clutter Return Map (dB)');
    % % % % % colorbar;
    % % % % % set(gca,'YDir','normal');
    % % % % % maxval = max(bindata(:));
    % % % % % caxis([maxval - 30 maxval]);
    % % % % % % bindata1 = bindata;
    
    
    
    
    % THIS WAY IS INDEED FASTER THAN METHOD 1, AND IT PRODUCES EXACT SAME
    % OUTPUT.
    % % disp('2st method...')
    bindata = zeros(length(rbin),length(dbin)) * nan;
    tic
    kr = discretize(clutter.rg, rbin);
    kd = discretize(clutter.dop, dbin);
    numFound = zeros(length(rbin),length(dbin));
    for i = 1:length(rbin)-1
        for j = 1:length(dbin)-1
            res = find( (kr == i) & (kd == j) );
            if(~isempty(res))
                numFound(i,j) = length(res);
                bindata(i,j) = max(clutter.power(res));
            end
        end
    end
    toc
    figure; add_print_callbacks;
    imagesc(1e-3*dbin,1e-3*rbin,bindata);
    hold on;
    C = contourc(1e-3*dbin,1e-3*rbin,numFound,[1 1]);
        % plot individual contours and color them WHITE
        %    cannot use contour() because it will map line color to existing
        %    colormap...
        k1 = 0;
        k2 = 0;
        ijk = 0;
        while k1 < size(C,2)
            ijk = ijk + 1;
            k1 = k1 + 1;
            k2 = C(2,k1);
            C(:,k1)
            [k1 k2]
            k1 = k1 + 1;
    
            x = C(1,k1+[0:k2-1]);
            y = C(2,k1+[0:k2-1]);
            plot(x,y,'w')
            k1 = k1 + k2 - 1;
        end
    
    % set(hcontour,'LineColor,','k')
    xlabel('Unambiguous Doppler (kHz)');
    ylabel('Unambiguous Range (km)');
    title('Clutter Return Map (dB)');
    colorbar;
    set(gca,'YDir','normal');
    maxval = max(bindata(:));
    caxis([maxval - 30 maxval]);
    % bindata2 = bindata;
    
    
    
    
    
    
    % bindata1(find(isnan(bindata1))) = 0;
    % bindata2(find(isnan(bindata2))) = 0;
    % diffbindata = bindata1-bindata2;
    % clc;
    % unequal_bindata = find(diffbindata ~= 0)
    % figure; imagesc(diffbindata); colorbar;
    
    
    
    % can I avoid the loop entirely???
    % disp('3rd method...')
    % bindata = zeros(length(rbin),length(dbin)) * nan;
    % kr = discretize(clutter.rg, rbin);
    % kd = discretize(clutter.dop, dbin);
    % tic
    
    
    
    
    % return
    
    
    %% Fold based on PRF
    c = 3e8;
    % numPulses = 600;
    
    cpi = 4e-3;
    prf = 550e3;
    prf = 95e3;
    pri = 1/prf;
    numPulses = floor(cpi/pri);
    
    % cpi = pri * numPulses;
    
    % duty = 0.33;
    % pw = duty * pri;
    
    pw = 0.5e-6;
    duty = pw/pri;
    
    rbin = c*(0:pw:(pri-pw))/2;
    dbin = 0:(prf/numPulses):prf;
    
    clutter.rg_amb = mod(clutter.rg, c*pri/2);
    clutter.dop_amb = sign(clutter.dop).*mod(abs(clutter.dop),prf);
    
    figure('Position',[72 162 1502 560]); add_print_callbacks;
    alpha_map = 10.^(clutter.power/10) ./ max(max(10.^(clutter.power/10)));
    subplot(1,2,1); him(1)=imagesc(1e-3*xg(1,:),1e-3*yg(:,1),clutter.rg_amb); title('Amb Rng @ Gnd')
    subplot(1,2,2); him(2)=imagesc(1e-3*xg(1,:),1e-3*yg(:,1),clutter.dop_amb); title('Amb Dop @ Gnd')
    % alpha(him(1),alpha_map);
    % alpha(him(2),alpha_map);
    
    bindata = zeros(length(rbin),length(dbin)) * nan;
    % % % tic
    % % % for kr = 1 : length(rbin) - 1
    % % %     for kd = 1 : length(dbin) - 1
    % % %         res = find( ...
    % % %             clutter.rg_amb > rbin(kr) & ...
    % % %             clutter.rg_amb < rbin(kr+1) & ...
    % % %             clutter.dop_amb > dbin(kd) & ...
    % % %             clutter.dop_amb < dbin(kd+1));
    % % %         if(~isempty(res))
    % % %             % get maximum for now
    % % %             bindata(kr,kd) = max(clutter.power(res));
    % % %         end
    % % %     end
    % % %     disp(kr./length(rbin));
    % % % end
    % % % toc
    tic
    kr = discretize(clutter.rg_amb, rbin);
    kd = discretize(clutter.dop_amb, dbin);
    numFound = zeros(length(rbin),length(dbin));
    for i = 1:length(rbin)-1
        for j = 1:length(dbin)-1
            res = find( (kr == i) & (kd == j) );
            if(~isempty(res))
                numFound(i,j) = length(res);
                bindata(i,j) = max(clutter.power(res));
            end
        end
    end
    toc
    figure; add_print_callbacks;
    imagesc(1e-3*dbin,1e-3*rbin,bindata);
    xlabel('Ambiguous Doppler (kHz)');
    ylabel('Ambiguous Range (km)');
    
    fmtstr = '%10.2f';
    title({'Clutter Return Map (dB)';...
        [num2str(prf/1e3,fmtstr),'kHz PRF | ',num2str(numPulses),' pulses | ',num2str(duty*100,fmtstr),'% duty | ',num2str(cpi*1e3,fmtstr),' ms dwell']});
    colorbar;
    set(gca,'YDir','normal');
    maxval = max(bindata(:));
    caxis([maxval - 30 maxval]);
    
    
    %% main beam data
    
    idx = find( clutter.antgain >= (clutter.antgain_peak - 25) ); % 3dB beamwidth
    
    % [idx1,idx2] = find( clutter.antgain >= (clutter.antgain_peak - 3) ); % 3dB beamwidth
    
    mbClutter.Dop    = clutter.dop(idx);
    mbClutter.DopAmb = clutter.dop_amb(idx);
    mbClutter.Rg     = clutter.rg(idx);
    mbClutter.RgAmb  = clutter.rg_amb(idx);
    
    % as1 = clutter.dop_amb(idx1,idx2);
    % as2 = clutter.rg_amb(idx1,idx2);
    
    figure('Position',[767 164 789 486]); add_print_callbacks;
    subplot(2,2,1); histogram(mbClutter.Dop*1e-3); title('MBC Dop'); xlabel('Dop [kHz]')
    subplot(2,2,2); histogram(mbClutter.DopAmb*1e-3); title('MBC Dop (Amb)'); xlabel('Dop [kHz]')
    subplot(2,2,3); histogram(mbClutter.Rg*1e-3); title('MBC Rng'); xlabel('Rng [km]')
    subplot(2,2,4); histogram(mbClutter.RgAmb*1e-3); title('MBC Rng (Amb)'); xlabel('Rng [km]')
    
    
    % figure;
    % histogram2(mbClutter.DopAmb,mbClutter.RgAmb,100)
    
    
    
    
    

