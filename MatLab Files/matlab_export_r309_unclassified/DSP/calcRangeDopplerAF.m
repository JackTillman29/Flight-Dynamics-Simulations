function varargout = calcRangeDopplerAF(sig,kernel,dt,varargin)
% Calculate the ambiguity function (range vs doppler) of the input signal
%
%   calcRangeDopplerAF(x,kernel,dt)
%   calcRangeDopplerAF(x,[],dt)
%       ** x is also the kernel
%    --> will generate a couple figures showing the AF
%       N = length(sig)
%       AF will be NxN
% [AF,Range,Doppler] = calcRangeDopplerAF(x,kernel,dt)
%
% [AF,Range,Doppler] = calcRangeDopplerAF(x,kernel,dt,dopExtent,nDop)
%       optional inputs (must be together):
%           dopExtent: bounds on doppler frequencies
%           nDop:      # of doppler 
%
% [AF,Range,Dopper] calcRangeDopplerAF(x,kernel,dt)
%    --> AF (NxN) matrix where N = length(x)
%    --> Range and Doppler (lxN)
% calcRangeDopplerAF(x,dt)
%    --> will generate a couple figures showing the AF
%        N = length(sig)
%        AF will be NxN
% [AF,Range,Dopper] = calcRangeDopplerAF(x,dt)
%    --> AF (NxN) matrix where N = length(x)
%    --> Range and Doppler (lxN)
%
% author: Jeff Hole (Booz Allen Hamilton), 2018

% ENSURE SIGNALS ARE IN ROW VECTORS
if(size(sig,1)>1)
    sig = sig.';
end
if(isempty(kernel) )
    kernel = sig;
else
    if(size(kernel,1) > 1)
        kernel = kernel.';
    end
end

% zero-pad signal and kernel with total length (need to compute FULL CONV)
N1 = length(sig);
N2 = length(kernel);
N2a = floor(N2/2);
N2b = N2 - N2a - 1;
% ZERO-PAD x WITH LENGTH-y ZEROS
% sig = [zeros(1,N2a) sig zeros(1,N2b)];
sig = [zeros(1,N2-1) sig];

c = 3e8;
% N = length(sig);
N = N1 + N2 - 1;
Fs = 1/dt;
if(length(varargin) == 0)
    dopArray = linspace(-Fs/2,Fs/2,N);
    nDop = N;
    M = N;
elseif(length(varargin) == 2)
    dopExtent = varargin{1};
    nDop = varargin{2};
    M    = nDop;
    if(length(dopExtent) == 1)
        dopExtent = [-1 1] * dopExtent;
    end
    dopArray = linspace(dopExtent(1),dopExtent(2),nDop);
else
    error('need to include both optional inputs')
end






t = [0:(N-1)]*dt;

K = length(kernel) ;
% % % sig and kernel are ROW VECTORS AT THIS POINT
% % af = zeros(N) ;
% % k = 0;
% % for dop = dopArray
% %     k = k + 1;
% %     disp(['k = ',num2str(k) ,' of ',num2str(nDop) ])
% %     tmp = PulseCompressionXCORRsStruct(sig .* exp(1i*pi*dop.*t) , ...
% %     kernel,'align') ;
% %     % AF first dim: Range/TimeDelay, second dim: Doppler
% %     % af (:, k) = tmp.signalCrossCor.';
% %     af(:,k) = tmp.signalCrossCorTrim(1:N).';
% % end

disp('Generating matrices...')
disp(['   size NxM = ',num2str(N),'x',num2str(M)])
if(N*M >= 8000^2)
    error(['this will take a long time!' ...
        ' go into code and change upper limit if you dare...'])
end

sigmat = repmat(sig.',[1 M]);

if(length(varargin) == 2)
    iarr = 1:N;
%     narr = dopArray/Fs * M;
    % to convert freq to N in cos(2*pi*i/N) == cos(2*pi*f*t)
    narr = Fs ./ dopArray;
else
%     iarr = ([0:(M-1)]-floor(M/2)) / M;
%     narr = 1:N;
    
    iarr = 1:N;
%     narr = ([0:(M-1)]-floor(M/2));
    narr = [M:-1:1] - floor(M/2);
%     narr = [M:-1:1];
end

[n,i] = meshgrid(narr,iarr);
% [i,n] = meshgrid([0:(M-1)],[1:N]);
%   i should be increasing along dim-2; looking at a single row, vals
%         should increase
%      i == frequency
%   n should be increasing along dim-1.
%      n == signal samples (time)
dopmat = exp(-1i*2*pi*i./n);
%     figure; imagesc(real(dopmat)); colorbar; title('inside calcRangeDopplerAF.m')
    % % figure;
    % % subplot(1,2,1); imagesc(i); colorbar
    % % subplot(1,2,2); imagesc(n); colorbar


%    dopmat = [s_fd0(1) s_fd1(1) ... s_fdN-1(1)]
%             [s fd0(2) s fd1(2) ... s fdN-1(2)]
%                 .        .            .
%                 .        .            .
%             [s_fdO(N) s_fd1(N) ... s_fdN-1(N)]
sigdopmat = sigmat .* dopmat;

kernelmat    = repmat(reshape(kernel,K,1),[1 M]);
% fftSigDopMat = fftshift(fft(sigdopmat,N,1),1); % convert each signal (column) into freq domain
% fftKernel    = fftshift(fft(kernelmat,N,1),1);
fftSigDopMat = fft(sigdopmat,N,1); % convert each signal (column) into freq domain
fftKernel    = fft(kernelmat,N,1);
%     figure; imagesc(real(kernelmat))

af = ifft( fftSigDopMat .* conj(fftKernel),N,1 );
% af = ifftshift( af,1 );

td = linspace(-floor(N/2),floor(N/2),N) * dt;

% normalize AF
af = af./max(max(abs(af)));

if(nargout > 0)
    varargout{1} = af;
    varargout{2} = td*c;
    varargout{3} = dopArray;
else
    % find sidelobe level relative to peak@ OdB
    % (mean of all points -3dB down)
    iTrim  = [floor(N/2):(N+floor(N/2))];
    afTrim = af(iTrim);
    idxSL  = find(10*log10(abs(afTrim)) <= -3);
    [idxML, idyML] = find(10*log10(abs(af)) > -3);
    
    mean_sidelobe_level = 10*log10(mean(abs(afTrim(idxSL))));
    
    mainlobe_rng_res = abs(td(max(idxML))*c);
    mainlobe_dop_res = abs(dopArray(max(idyML)));
    
    figpos = [694 303 793 592];
    
    % RANGE VS DOPPLER FREQ
    figure('Position',figpos); add_print_callbacks;
    imagesc(c*td, dopArray*1e-3, 10*log10(abs(af.')));
    hcb = colorbar;
    set(get(hcb,'YLabel'),'String','dB','FontSize',16,'FontWeight','bold')
    caxis([-20 0])
    set(gca,'YDir','normal')
    xlabel('Range [m) ')
    ylabel('Doppler [kHz] ' )
    title({'Ambiguity Function'; ...
        ['\color{gray}\DeltaR = ',num2str(mainlobe_rng_res),'m', ...
        ' | \DeltaF = ',num2str(mainlobe_dop_res*1e-3), ...
        'kHz', '\color{black}']; ...
        ['\color{gray}Avg sidelobe level: ', ...
            num2str(mean_sidelobe_level),'dB\color{black}']})

end

end
