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
%aesa.aVar = 0*0.5;
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



%% create a noise source out there
nSamples = 1000;
nSubArrays = 16;
close all;

nSteps = 1000; % jam steps
nAngSamples = 5000;
angSpan = [-30 10]*pi/180;%-30, 10
outArr = nan*zeros(nSteps,nAngSamples);
jam_az = linspace(angSpan(1)*180/pi,angSpan(2)*180/pi,nSteps);
angleSteerDeg = [-20 0];
for kjam = 1 : length(jam_az)
    jam_az(kjam)
    noise_azimuths = [jam_az(kjam)]*pi/180;
    noise_elevations = [0]*pi/180;
    
    
    
    % Apply full array phase shifting
    aesa.upd_setSteer(angleSteerDeg(1)*pi/180,angleSteerDeg(2)*pi/180);
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
    
    % Now perform MUSIC processing
    % Compute MUSIC for array of far-field rel main beam
    nAz = 400;
    nEl = 400;
    azVec = linspace(-10,10,nAz);
    elVec = linspace(-10,10,nEl);
    
    % sh (sha) will be [N subarrays x (N Az x N El)]
    sha = conj( exp(1i*aesa.upd_setSteerArray( ...
                azVec*pi/180, ...
                elVec*pi/180, ...
                'subarray')));
            
    % project all steering vectors onto noise subspace (method 2)
    % NOTE: THis is done slightly differently than above because x'*x does not
    % yeild a scalar dot product for the full matrix implementation. This
    % requires the dot product to be implemnted as an element-wise multiply and
    % summation. They are equivalent.
    tmp = (eye(nSubArrays)-V_is*V_is')*sha;
    
    noiseSubProjArr2 = reshape(sum(conj(sha).*tmp,1),nEl,nAz);
    
    MUSIC_Arr = 10*log10(1./abs(noiseSubProjArr2));
    if(kjam == 1)
        figure;
        him=imagesc(azVec,elVec,MUSIC_Arr);
        xlabel('\DeltaAzimuth Angle of Arrival');
        ylabel('\DeltaElevation Angle of Arrival');
        title('MUSIC Az/El Response [dB]');
        colorbar;
        set(gca,'YDir','normal');
        set(gcf,'Position',[1137          48         560         420]);
        %caxis([-10 10]);
    else
        set(him,'CData',MUSIC_Arr);
    end
    
%     return;
end