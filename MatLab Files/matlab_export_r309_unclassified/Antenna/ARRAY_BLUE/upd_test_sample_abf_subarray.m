close all;
clear all;
clc;

% Notes:
%  N = # of elements

% Instantiate array object
aesa = array('Sample Array');

% Initialize the array object
aesa.InitStapArray( ...
    3e8/34e9, ...               % wavelength, in meters
    1.0, ...                    % efficiency, between 0.0 and 1.0
    'Sample_Elements.txt', ...  % element position file, x,y,z meters [3 x N]
    [], ...                     % amplitude sum taper file, [N x 1]
    [], ...                     % amplitude delta az taper file, [N x 1]
    [], ...                     % amplitude delta el taper file, [N x 1]
    []);                        % element pattern file, not implemented in InitStapArray

% Create the array taper using internal routines
aesa.compTaylorWgt(4,-30);      % compute taylor (rx) weights
aesa.compBaylissWgt(4,-40);     % compute difference pattern weights

% set amplitude & phase noise
aesa.aVar = 0*0.5;
aesa.pVar = 1*10*pi/180;

% set element pattern (NOTE: This simply creates a 1-D solid angle cos^2
% element pattern with a peak power gain of elGainDbi).
elGainDbi = 4.5;
elGainDbi_V = elGainDbi / 2.0;
elGain_V = 10^(elGainDbi_V/10);
aesa.matlab_setElementPattern( ...
    linspace(-pi/2,pi/2,100), ...
    elGain_V .* cos(linspace(-pi/2,pi/2,100)));


aesa.upd_redefineCoordinates();
aesa.upd_scaleWeights();


%% Create sub-array definition
%aesa.upd_subArrayNumber = ...
%    [ones(1,32 )*1 ones(1,96)*2 ones(1,96)*3 ones(1,128)*4 ...
%     ones(1,128)*5 ones(1,96)*6 ones(1,96)*7 ones(1,32 )*8 ]';

% sub array
load new_sub_original.mat;
for ksub = 1 : 16
    %    yd = eval(sprintf('new_sub%d(:,2);',ksub));
    subArrayIdx = eval(sprintf('new_sub%d(:,1);',ksub));
    % hp=plot(yd,zd,'.','MarkerSize',18);
    % disp(length(zd))
    aesa.upd_subArrayNumber(subArrayIdx) = ksub;
end
clear subArrayIdx new_sub*;

aesa.upd_setSubArrayWeight(ones(16,1));
aesa.plotGeometry(aesa.upd_subArrayNumber);
return;


%% create a noise source out there
nSamples = 1000;
nSubArrays = 16;
close all;
%aesa.pVar = 5*pi/180;
nSteps = 1000;
nAngSamples = 5000;
angSpan = [-50 10]*pi/180;
outArr = nan*zeros(nSteps,nAngSamples);
jam_az = linspace(angSpan(1)*180/pi,angSpan(2)*180/pi,nSteps);
for kjam = 1 : length(jam_az)
    
    noise_azimuths = [jam_az(kjam)]*pi/180;
    noise_elevations = [0]*pi/180;
    
    
    
    % Apply full array phase shifting
    aesa.upd_setSteer(-20*pi/180,0*pi/180);
    %if(kjam == 1)
    subArray_iq = 0;
    for k = 1 : length(noise_azimuths)
        noise_sig = randn(1,nSamples) + 1i * randn(1,nSamples);
        %noise_sig = 1 + 1i*0;
        noise_azimuth = noise_azimuths(k);
        noise_elevation = noise_elevations(k);
        
        % clear sub-array weights
        aesa.upd_setSubArrayWeight(1);
        
        subArray_iq = subArray_iq + ...
            aesa.upd_genIQ_Channels( ...
            noise_sig, ...
            'fromAz', noise_azimuth, ...
            'fromEl', noise_elevation, ...
            'subArray', 1:nSubArrays, ...
            'applySteering', 1, ...
            'sumOutput', 0);
    end
    %end
    
    % Estimate correlation matrix
    Rest = aesa.upd_estCorrMtx(subArray_iq,[]);
    
    %
    [Vc,Dc]=eig(Rest);
    [~,sortIdx]=sort(Dc,'descend');
    Dc = diag(Dc(sortIdx));
    Vc = Vc(:,sortIdx);
    eigIdx = length(noise_azimuths);
    V_is = Vc(:,1:eigIdx);
    
    % *** what is the s vector for the sub array??
    %     1 - is it just ones vector (steering is already done up in the
    %     array yes!
    %     2 - does it include the steering vector? NO!
    nSubArrays = length(unique(aesa.upd_subArrayNumber));
    sh = ones(nSubArrays,1);
    shprime = (eye(nSubArrays)-V_is*V_is') * sh;
    abf_termc = conj(shprime);
    
    
    
    if(kjam == 1)
        aesa.upd_setSubArrayWeight(ones(nSubArrays,1));
        hf = figure(1);
        hp = aesa.upd_plot(linspace(angSpan(1),angSpan(2),nAngSamples),0,aesa.ENUM_MODE_UNIFORM);
        set(hp,'Color',0.4*[1 1 1]);
        ha = gca;
        hL = line([0 0],[-40 40],'Color','r');
    end
    
    aesa.upd_setSubArrayWeight(abf_termc);
    
    if(kjam == 1)
        hold on;
        [hp]=aesa.upd_plot(linspace(angSpan(1),angSpan(2),nAngSamples),0,aesa.ENUM_MODE_UNIFORM);
        set(hp,'LineWidth',2);
        plotData = get(hp,'YData');
        set(gcf,'Position',[6    47   560   420]);
        figure;
        hi=imagesc(angSpan*180/pi,jam_az,outArr);
        colormap(jet)
        caxis([-10 33]);
        xlabel('Yaw Plane [deg]');
        ylabel('Yaw EA AoA [deg]');
        title('Adaptive Beamformer Response');
        colorbar;
        set(gcf,'Position',[571    47   560   420]);
        
    else
        [~,plotData]=aesa.upd_plot(linspace(angSpan(1),angSpan(2),nAngSamples),0,aesa.ENUM_MODE_UNIFORM);
        set(hp,'YData',plotData);
    end
    set(hL,'XData',[0 0] + jam_az(kjam));
    drawnow;
    
    outArr(kjam,:) = plotData;
    set(hi,'CData',outArr);
    
%     %% temporarily undo abf
%     close all;
%     aesa.upd_setSubArrayWeight(1);
%     
%     aesa.upd_computeSubArrayCenters();
%     
%     saPos = aesa.upd_subArrayPos;
%     
%     aesa.upd_setSteer(0*pi/180,0*pi/180);
%     
%     % Build steering vector components
%     azSteer = 0;
%     elSteer = 0;
%     
%     [uxs,uys,uzs]=sph2cart( ...
%         azSteer * pi/180, ...
%         elSteer * pi/180, ...
%         1 ); % az,el,1
%     
%     % Form steering vector
%     rs = [uxs uys uzs]';
%     
%     % Compute phasing (rad) over array elements [N x 1]
%     sa_phs_steer = -aesa.twopi_ovr_lambda * saPos * rs;
%     
%     aesa.upd_setSubArrayWeight(exp(1i*sa_phs_steer).');
%     aesa.upd_plot();
%     
%     %%
%     testAngleAz =  2 * pi/180;
%     testAngleEl =  0 * pi/180;
%     
%     nSigEig   = 1;
%     nNoiseEig = length(Vc(:,1)) - nSigEig;
%     
%     [uxs,uys,uzs] = sph2cart( ...
%         testAngleAz * pi/180, ...
%         testAngleEl * pi/180, ...
%         1 ); % az,el,1
%     
%     rs_music = [uxs uys uzs]';
%     
%     sa_music_phs_steer = -aesa.twopi_ovr_lambda * saPos * rs_music;
%     sa_music_wgt = exp(1i*sa_music_phs_steer);
%     
%     sa_music_wgt.'*(eye(16)-Vc(:,1)*Vc(:,1)')'*sa_music_wgt
%     
    
%     return;
end