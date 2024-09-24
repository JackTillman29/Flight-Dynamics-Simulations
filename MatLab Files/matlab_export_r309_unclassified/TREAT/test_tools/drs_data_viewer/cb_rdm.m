function cb_rdm( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
uiwait(msgbox('Click the start and stop points','Region'));
rdm_window = ginput(2);
ezja_header = evalin('base','ezja_header;');

uiwait(msgbox('Click the threshold for reference pulse present','Level'));
level = ginput(1);
level = level(2);

    haxes = findobj('Type','Axes','Parent',gcf);
    ref_channel = evalin('base','ref_channel;');
    
    hax = haxes(1);
    xb=get(hax,'XLim');
    xvec = evalin('base','xvec;');
    dat1 = evalin('base','dat1;');
    dat2 = evalin('base','dat2;');
    
    %iStart = find(xvec > xb(1),1,'first');
    %iStop = find(xvec < xb(2),1,'last');
    
    iStart = find(xvec > rdm_window(1,1),1,'first');
    iStop = find(xvec < rdm_window(2,1),1,'last');
    
    if(~ezja_header)
        i1_counts = dat1.pageMemory(iStart:iStop,1);
        q1_counts = dat1.pageMemory(iStart:iStop,2);
        i2_counts = dat2.pageMemory(iStart:iStop,1);
        q2_counts = dat2.pageMemory(iStart:iStop,2);

        v1 = single(i1_counts) + 1i*single(q1_counts);
        v2 = single(i2_counts) + 1i*single(q2_counts);
    else
         dat1.temp = single(dat1.pageMemory(iStart:iStop));
         dat2.temp = single(dat1.pageMemory(iStart:iStop));
         disp('Applying Hilbert Transform on Channel 1,2');
         dat1.tempf = fft(dat1.temp);
         dat2.tempf = fft(dat2.temp);
         dat1.tempf(floor(length(dat1.tempf)/2):end) = 0;
         dat2.tempf(floor(length(dat2.tempf)/2):end) = 0;
         v1 = ifft(dat1.tempf);
         v2 = ifft(dat2.tempf);
         disp('Done');
        
    end
    
    f1 = fftshift(fft(v1));
    f2 = fftshift(fft(v2));
    
    t = evalin('base','Ts;')* (1:length(v1))';
    
    fvec = evalin('base','Fs;') * (0:length(f1)-1)./length(f1);
    

    
    if(ref_channel == 1)
        [fmax,fmaxi]=max(abs(fft(v1)));
        disp('Using channel 1 for reference');
    else
        [fmax,fmaxi]=max(abs(fft(v2)));
        disp('Using channel 2 for reference');
    end
    
    coho_est_value = fvec(fmaxi);
    disp(['Estimated Center Frequency ' num2str(coho_est_value)]);
    cohoData = exp(1i*2*pi*(coho_est_value)*t)';
    
    
    ref_channel_dop_offset = evalin('base','ref_channel_dop_offset;');
    cohoDataMod = exp(1i*2*pi*(ref_channel_dop_offset)*t)';
    if(ref_channel_dop_offset ~= 0)
        warning(['Using ref channel offset of ' num2str(ref_channel_dop_offset) ' Hz']);
    end
    
    %%
    inData.dt = evalin('base','Ts;');
    inData.t = xvec;
    
    outData.dt = evalin('base','Ts;');
    outData.t = xvec;
    
    if(ref_channel == 1)
        inData.Vc = v1 ./ cohoDataMod(1:length(v1)).';
        outData.Vc = v2;
    else
        inData.Vc = v2 ./ cohoDataMod(1:length(v1)).';
        outData.Vc = v1;
    end
    
    
%%
inSigAboveThresh = 20*log10(abs(inData.Vc)) >  level;
idxPulseStart = find(diff(inSigAboveThresh) == 1);
idxPulseStop  = find(diff(inSigAboveThresh) == -1);

blank = 0;
if(blank)
    blank_vec = inSigAboveThresh;
end

% fix trim issues
% requirements: start(1) < stop(1), and stop(end) > start(end)
if(idxPulseStart(1) > idxPulseStop(1))
    idxPulseStop = idxPulseStop(2:end);
end
if(idxPulseStart(end) > idxPulseStop(end))
    idxPulseStart = idxPulseStart(1:end-1);
end


% Quick mod to center
% get half a pulsewidth
nHalfPulse = idxPulseStop(1)-idxPulseStart(1);
%if(nHalfPulse < 0) % caught the end of a pulse in the beginning
%    idxPulseStop = idxPulseStop(2:end);
%end


if(nHalfPulse > 0)
    idxPulseStart = idxPulseStart - nHalfPulse;
    idxPulseStop = idxPulseStop - nHalfPulse;
    
    if(idxPulseStart(1) < 1)
        idxPulseStart = idxPulseStart(2:end);
        idxPulseStop = idxPulseStop(2:end);
    end
    
end


nPulseStart = length(idxPulseStart);
fprintf('Estimated # of Pulses: %d\n',nPulseStart);

priPerLine = mean(diff(idxPulseStart)) * inData.dt;
prfEst = 1./priPerLine; % for later
fprintf('Estimated PRI (us): %18.6f\n',1e6*priPerLine);
priPerLine = 1.05 * priPerLine;
smpPerLine = floor(priPerLine ./ inData.dt);

prf_str = sprintf('%6.2fkHz (%6.2f\\mus) ',1e-3*prfEst,1e6/prfEst);

%
% repackage data
Min = complex(zeros(nPulseStart,smpPerLine,'single'));
Mout = complex(zeros(nPulseStart,smpPerLine,'single'));
Mcoho = complex(zeros(nPulseStart,smpPerLine,'single'));
Mpp = zeros(nPulseStart,smpPerLine,'single');
fastTimeVec = (0:size(Min,2)-1) * inData.dt * 1e6;
slowTimeVec = (inData.t(idxPulseStart) - inData.t(1)) * 1e-3;


for k = 1 : nPulseStart
    try
        Min(k,:) = inData.Vc(idxPulseStart(k) : idxPulseStart(k) - 1 + smpPerLine);
        Mout(k,:) = outData.Vc(idxPulseStart(k) : idxPulseStart(k) - 1 + smpPerLine);
        Mcoho(k,:) = cohoData(idxPulseStart(k) : idxPulseStart(k) - 1 + smpPerLine);
    catch
        if(k == nPulseStart)
            n = length(inData.Vc(idxPulseStart(k) : end));
            Min(k,1:n) = inData.Vc(idxPulseStart(k) : end);
            Mout(k,1:n) = outData.Vc(idxPulseStart(k) : end);
            Mcoho(k,1:n) = cohoData(idxPulseStart(k) : end);
        else
            warning(['Data remap error on pulse ' num2str(k) ' of ' num2str(nPulseStart)]);
        end
    end
end

Min = Min .* Mcoho;
Mout = Mout .* Mcoho;

nColumns = size(Min,2);
hColumns = floor(nColumns*0.25);

%Min = Min(:,[hColumns:


figure;
subplot(1,2,1);

imagesc(fastTimeVec,slowTimeVec,20*log10(abs(Min)));
xlabel('PRI Time (\mus)');
ylabel('CPI Time (ms)');
title('Reference - FTST');
ah = gca;
colorbar;
subplot(1,2,2);
imagesc(fastTimeVec,slowTimeVec,20*log10(abs(Mout)));
xlabel('PRI Time (\mus)');
ylabel('CPI Time (ms)');
title('Response - FTST');
ah = [ah gca];
colorbar;
%colormap(colormap_fade2black(256,.5));
linkaxes(ah);

set(gcf,'Position',[737          71        1169         443]);
PrepForPrint(gcf,evalin('base','pp;'));
add_print_callbacks;

figure;
subplot(1,2,1);
dvec = linspace(-prfEst/2,prfEst/2,length(slowTimeVec));
rdm_1 = 20*log10(abs(fftshift(fft(Min,[],1))).');
rdm_2 = 20*log10(abs(fftshift(fft(Mout,[],1))).');
max_1 = max(rdm_1(:));
max_2 = max(rdm_2(:));
max_both = max([max_1 max_2]);

imagesc(dvec*1e-3,fastTimeVec,rdm_1);
ylabel('PRI Time (\mus)');
xlabel('Relative Doppler (kHz)');
title('Reference - RDM');
ah = gca;
%caxis([max_both-60 max_both]);
colorbar;

subplot(1,2,2);

imagesc(dvec*1e-3,fastTimeVec,rdm_2);
ylabel('PRI Time (\mus)');
xlabel('Relative Doppler (kHz)');
title('Response - RDM');
ah = [ah gca];
%caxis([max_both-60 max_both]);
colorbar;
linkaxes(ah);

colormap(colormap_fade2black(256,.5));
set(gcf,'Position',[736         598        1171         379]);
PrepForPrint(gcf,evalin('base','pp;'));
add_print_callbacks;
if(0)
    fosvec = linspace(-0.5,0.5,size(Mout,2)*10);
    figure;
    imagesc(1e-3*fosvec*evalin('base','Fs;'),slowTimeVec,20*log10(abs(fftshift(fft(Mout,10*size(Mout,2),2)))));
    title('Waterfall');
    xlabel('kHz');
    ylabel('ms');
end
end

