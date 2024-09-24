function [...
    inputCrossCorr,...
    inputCrossCorrTrim,...
    kleadlag,...
    kernelAutoCorr,...
    fftInput,fftKernel,fftFreqVect] = PulseCompressionXCORRs(input,kernel,timealign,varargin)
% varargin{1} = output unit scale factor of x axis. If omitted, it is considered
% "samples
debug = 0;

if(size(input,2) > size(input,1))
    row_vec_data = 0;
    input = input.';
    if(size(kernel,2) > size(kernel,1))
        kernel = kernel.';
    end
else
    row_vec_data = 1;
    if(size(kernel,2) > size(kernel,1))
        kernel = kernel.';
    end
end




% code expects column vectors...
input = input.';
kernel = kernel.';

% NOTES; Appears that the "shift" needed to make the data look correct is
%        based on the kernel length only.

% K. Sawmiller 2012
    if(~exist('timealign'))
        timealign = 'align';
        %disp('Defaulting to time-aligned signal (input will be padded to be atleast size of kernel)');
    end
    
    if(strcmp(timealign,'align'))
        padlen = round( (length(kernel)-length(input))/2 );
        if(padlen > 0)
            input = [0*(1:padlen) input 0*(1:padlen)];
        end
    end

    % Lengths and sizes. The length of the convolved signal will be
    N = length(input );
    M = length(kernel);
    L = N + M - 1;
    if(debug)
        fprintf('N=%d, M=%d, L=%d\n',N,M,L);
    end

    % Subtract any residual mean (DC) off the signals
    % ** This step has been removed. If needed, it should be applied to the
    % input **
    
    % Compute the FFT of each waveform
    fftInput  = fft(input  ,L);
    fftKernel = fft(kernel ,L);
    
    % Fourier theory says that the FFT of the cross correlation of two signals in the
    % time domain is equal to the product of one signal FFT times the
    % complex conjugate of the other signal FFT
    
    % Determine correlation values for each signal
    kernelAutoCorr = ifft(  fftKernel .* conj( fftKernel  ) );
    inputCrossCorr = ifft(  fftInput  .* conj( fftKernel  ) );
    kernelAutoCorr = ifftshift(kernelAutoCorr)./(sqrt(M)*sqrt(M));
 
%     scalefac = L * M;
%     temp_Xm = inputCrossCorr ./ scalefac;
%     Xm(1:floor(M/2)) = temp_Xm((L-floor(M/2)+1):L);
%     Xm((floor(M/2)+1):N) = temp_Xm(1:(N-floor(M/2)));
%     inputCrossCorr = Xm;
    % After inputing a pulse train, it appears that there is a full copy of the
    % kernel length that needs to be moved in front (fft shift does not do
    % this correctly)
     if(1)
         % first, "unwrap" the data (full kernel worth of data time reversed on end)
         % Reference: M = Kernel length, N = Sample Length, L = Convolution Length (M+N-1)
         
         %% AMIR HACKS
         inputCrossCorr = [ ...
             inputCrossCorr((L-M+1):end) ...
             inputCrossCorr(1:(L-M)) ];
         
         % next, "trim" the data. Half the kernel needs trimmed on both sides
         K2 = floor(M/2);
         
         %inputCrossCorrTrim = inputCrossCorr(K2:(end-K2-1));
         
         inputCrossCorrTrim = inputCrossCorr(K2:(end-K2));
         
         % because of the flooring function, sometimes the trimmed length
         % can be off by one sample from the input
         if(length(inputCrossCorrTrim) ~= length(input))
             inputCrossCorrTrim = inputCrossCorr(K2:(end-K2-1));
         end
         
         % scale the data based on sample length
         inputCrossCorrTrim = inputCrossCorrTrim ./ M;
         inputCrossCorr     = inputCrossCorr ./ M;
         
                  %% AMIR HACKS
%          inputCrossCorrTrim =inputCrossCorr(1:L);  ALK commented out??
     %inputCrossCorr = [ ...
     %    inputCrossCorr((end-floor(M)):end) ...
     %    inputCrossCorr(1:(end-floor(M)-1)) ...
     %    ] ./ ...
     %    (sqrt(M)*sqrt(N)) * sqrt(N/M);
     end
     
    %inputCrossCorr = fftshift(inputCrossCorr)./(sqrt(M)*sqrt(N)) * sqrt(N/M);
    
    fftFreqVect    = (0:(L-1))./L;
    
    kleadlag = (-L/2):(L/2);
    if( length(kleadlag) > L )
        kleadlag = kleadlag(2:end);
    else
        kleadlag = [kleadlag kleadlag(end)+1];
    end
    
    % now apply x axis units if provided
    if(nargin == 4 && ~isempty(varargin{1}))
        kleadlag = kleadlag .* varargin{1}{1};
    end
    
    if(nargout == 0)
        figure;
        subplot(2,1,1);
        plot( ...
            kleadlag,real(kernelAutoCorr),'b:', ...
            kleadlag,abs(kernelAutoCorr),'r', ...
            'LineWidth',2);
        title({'\fontsize{12}\bf{}Reference Signal Auto-Correlation'});
        xlabel({'\bf{}Sample'});
        ylabel({'\bf{}Magnitude'});
        grid on;
        legend( ...
            '\fontsize{8}\bf{}True Signal', ...
            '\fontsize{8}\bf{}Envelope', ...
            'location','best');
        subplot(2,1,2);
        plot( ...
            kleadlag,real(inputCrossCorr),'b:', ...
            kleadlag,abs(inputCrossCorr),'r', ...
            'LineWidth',2);
        title({'\fontsize{12}\bf{}Input Signal Cross-Correlation'});
        xlabel({'\bf{}Sample'});
        ylabel({'\bf{}Magnitude'});
        grid on;
        legend( ...
            '\fontsize{8}\bf{}True Signal', ...
            '\fontsize{8}\bf{}Envelope', ...
            'location','best');
    end
    
    if(row_vec_data)
        inputCrossCorr = inputCrossCorr.';
        inputCrossCorrTrim = inputCrossCorrTrim.';
        kleadlag = kleadlag.';
        kernelAutoCorr = kernelAutoCorr.';
        fftInput = fftInput.';
        fftKernel = fftKernel.';
        fftFreqVect = fftFreqVect.';
    end
    
end