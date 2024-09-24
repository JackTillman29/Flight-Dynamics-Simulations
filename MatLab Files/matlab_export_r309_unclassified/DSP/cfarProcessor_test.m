close all; clear all; clc

if(isunix)
    matlab_root = '/home/holeja/';
elseif(ispc)
    matlab_root = 'I:\MATLAB\'
end

addpath([matlab_root 'Printing'])
addpath([matlab_root 'Analysis'])
addpath([matlab_root 'DSP'])
addpath([matlab_root 'Signal Generation\Noise Signal'])
c = lines;

pf = @(x) 20*log10(abs(x));

% assumes single pulse (no pulse integration) CA-CFAR
%  CA == Cell Averaging
%  CFAR == Constant False Alarm Rate
% 

nTrain = 10;
nGuard = 1;
Pfa    = 1e-5;  % probability of false alarm

BW = 200e3;
Fs = 10*BW; dt = 1/Fs;

N = 100*(2*nTrain);
itgt = floor(N/2);
t = [1:N];

snr_dB = 10;
    k = 1.38e-23; T = 290; Nf = 10^(0/10);
    noisePSD = k * T * Nf;

nrep = 5000;
numDets    = 0;
numTrueDet = 0;
numFA      = 0;
% Pdet = zeros(1,nrep);
% Pfa = zeros(1,nrep);
% for irep = 1:nrep
    x = sqrt(noisePSD * Fs * 0.5) * (randn(1,N) + 1i*randn(1,N));
    %     wx = WaveObj(x,[],Fs,length(x));
    %     plot(wx)
    x(itgt) = x(itgt) + sqrt(noisePSD * Fs * 10^(snr_dB/10));


    % process signal (square law)
    x2 = abs(x);
        nCells = 2*nTrain;
        alpha = nCells .* (Pfa^(-1/nCells) - 1);
        alpha = 10^(8/20);

    % run CFAR on x2
    [y,thr,det] = cfarProcessor(x2,alpha,nTrain,nGuard);
        % ignore dets on edges
        det(1:(nTrain+nGuard)) = 0;
        det((end-(nTrain+nGuard)):end) = 0;
        
    figure; add_print_callbacks;
    plot(pf(x2))
    hold on;
    plot(pf(thr))
    plot(pf(y),'k.','MarkerSize',12)
    legend('signal envelope','threshold','detections')
    xlabel('Sample Number')

    idet = find(det == 1);

    numDets    = numDets + sum(det);
    numTrueDet = numTrueDet + length( find(t(det==1) == itgt) );
    numFA      = numFA      + length( find(t(det==1) ~= itgt) );

%     Pdet(irep) = numTrueDet/irep;
%     Pfa(irep)  = numFA/(irep*N);


%     if(mod(irep,200)==0)
%     if(~exist('hfig'))
%         hfig = figure('Position',[643 88 932 704]); add_print_callbacks;
%         hax(1) = subplot(4,1,1:2);
%         hp(1) = plot(t,abs(x));
%         hold on;
%         hp(2) = plot(t(itgt),abs(x(itgt)),'r.');
%         hp(3) = plot(t,thr,'color',c(2,:));
%         ylims = ylim(hax(1));
%         hp(4) = plot(t,det*ylims(2),'color',c(3,:))
% 
%         hax(2) = subplot(4,1,3);
%         hp(5) = plot(1:irep,Pdet(1:irep)*100); title('Probability of Detection [estimate]')
%         hax(3) = subplot(4,1,4);
%         hp(6) = plot(1:irep,Pfa(1:irep)*100); title('Probability of False Alarm [estimate]')
%     else
%         set(hp(1),'YData',abs(x))
%         set(hp(2),'XData',t(idet),'YData',abs(x(idet)))
%         set(hp(3),'YData',thr)
%         ylims = ylim(hax(1));
%         set(hp(4),'YData',det*ylims(2))
% 
%         set(hp(5),'XData',1:irep,'YData',Pdet(1:irep))
%         set(hp(6),'XData',1:irep,'YData',Pfa(1:irep))
% 
%     end



%     % title(hax(1),{['# True Dets: ',num2str(numTrueDet), ' / ',num2str(irep),' (# missed: ',num2str(irep - numTrueDet),', P_{det}: ',num2str(Pdet*100,'%10.1f'),' %)'];...
%     %     ['# False Alarms: ',num2str(numFA),' | Total # Cells Checked: ',num2str(irep * N),' (P_{fa}: ',num2str(Pfa*100,'%10.3f'),'%)']})
% 
%     title(hax(1),{['# True Dets: ',num2str(numTrueDet), ' / ',num2str(irep),' (P_{det}: ',num2str(Pdet(irep)*100,'%10.1f'),' %)'];...
%         ['# False Alarms: ',num2str(numFA),' | (P_{fa}: ',num2str(Pfa(irep)*100,'%10.3f'),'%)']})
% 
%     pause(0.1)
% 
%     xlim(hax(1),itgt + 200*[-1 1])
%     end

% end




