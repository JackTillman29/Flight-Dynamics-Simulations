

function [x, nIts, rOut] = cgToep(AStick, b, x, type, tol, nIts, useNfft, useGpu)
%  Runs CG, Assuming Toeplitz Matrix Input
%
% AStick: (M + N - 1)x1, where M >= N, Outside Edge of A matrix
% b: Mx1 Number of Data Points
% x: Nx1 Initial Weight Vector
%
%  JAH:
%  x = inv(A^H * A) * A^H * b
%       AStick == A^H (H == hermitian = conj(A'))
%
% Type:     0 -> Toeplitz
%           1 -> Block Toeplitz
%               Can be passed in as structure or Matrix (even blocks are
%               evenly sized
% tol: Tolerance (Default 10^-12)
% nIts: Number of Iterations (Default: size(A, 2))

debugFlag = 0;

if ~exist('type', 'var') || isempty(type)
    type = 0; % Default to Non-Square
end

if ~exist('tol', 'var') || isempty(tol)
    tol = 10^-12; % Default to Tolerance
end

if ~exist('useNfft', 'var') || isempty(useNfft)
    useNfft = 0;
end

if ~exist('useGpu', 'var') || isempty(useGpu)
    useGpu = 0;
end

switch type
    case 0  % This is All Toeplitz
        if size(AStick, 2) ~= 1
            error('''AStick'' MUST be outside edge of ''A'' Matrix');
        end
        
        M = size(b, 1);         % # of rows in A matrix
        N = size(AStick, 1) - M + 1; % # of delay taps in A matrix
        nStick = M + N - 1;
        
        if ~exist('x', 'var') || isempty(x)
            x = zeros(N, 1); % Initialize to all zeros
        end
        
        if ~exist('nIts', 'var') || isempty(nIts)
            nIts = N; % Default to Max Rank
        end
        rOut = zeros(1, nIts + 1);
        
        %nfft = 2^(ceil(log2(M + N - 1))); % Number of Points in FFT
        nfft = M + N - 1;
       
        if useGpu == 1
            AStick = gpuArray(AStick);
            b = gpuArray(b);
            rOut = gpuArray(rOut);
            x = gpuArray(x);
        end
        
        AStickFC = conj(flipud(AStick)); % Conj and flip so convolution == xcorr (A')
        
        if useNfft == 1
            fftA = fft(AStick, nfft);     % (A)
            fftAFC = fft(AStickFC, nfft); % (A')
        else
            %fftA = fft(AStick);                       % (A)
            %fftAFC = fft(AStickFC);                   % (A')
            fftA = fft([AStick; zeros(nfft-nStick, 1)]);     % (A)
            fftAFC = fft([AStickFC; zeros(nfft-nStick, 1)]); % (A')
        end
        
        % % This Computes an XCORR
        %   or = b;
        %   otemp = ifft(fftAFC.*fft(or,nfft)); % xcorr A with b (A'*b)
        %   op = otemp(M+[0:N-1]);              % Grab taps we care about (A'*b)
        %   gamma1 = op'*op;                    % Compute Power b'AA'b
        %
        
        %------------------------------------------------------------------
        %
        %  % -------------------- CODE DOES THIS USING FFTs ----------------
        %
        %  r = b - A*x;   % Initial Residual (Stuff in A that x hasn't solved for yet) --> This is negative gradient (toward solution)
        %  p = A'*r;      % Use initial Residual as our initial update guess (inital basis vector, orthogonal to intital guess Xo)
        %  rsold = p'*p;  % Power of Residual Vector
        %------------------------------------------------------------------
        if useNfft == 1
            temp = ifft(fftA.*fft(x, nfft));            % A*x
        else
            %temp = ifft(fftA.*fft(x));            % A*x
            temp = ifft(fftA.*fft([x; zeros(nfft-N, 1)]));            % A*x
        end
        r = b - temp(N+[0:M-1]);                    % r = b - A*x , length M
        if useNfft == 1
            temp = ifft(fftAFC.*fft(r,nfft));           % xcorr A with r (A'*r)
        else
            %temp = ifft(fftAFC.*fft(r));           % xcorr A with r (A'*r)
            temp = ifft(fftAFC.*fft([r; zeros(nfft-M, 1)]));           % xcorr A with r (A'*r)
        end
        p = temp(M+[0:N-1]);                        % Grab taps we care about (p = A'*r)
        rsold =p'*p;
        
        rOut(1) = rsold;                            % Hold Residuals
        
        for ii=1:nIts                                   % Loop over all dimensions
            %-----------------------------------------------------------------------------------------------------------------------------------------------------------------
            %
            %  % -------------------- CODE DOES THIS USING FFTs ----------------
            %   %Ap = A*p;                 % Projection of A onto P (How much of A lives on P)
            %   %alpha = rsold/(Ap'*Ap);   % Ratio of Power in Residual / Power of P,  A Matrix Projection (scales P properly so Ax =b)
            %   %x = x + alpha*p;          % Update Solution (Add this iterations orthogonal guess/basis (P) to existing)
            %   %r = r - alpha*Ap;         % Project this guess out of residual --> This is negative gradient (towards solution)
            %   %rsnew = (r'*A)*(A'*r);    % Computed New Residual (trying to drive this to zero --> means Ax = b
            %   %p = A'*r + rsnew/rsold*p; % Contains all n orthogonal basis vectors (the Ps) (with there powers normalized to power of newest basis vector)
            %   %rsold = rsnew;            % update power number
            %-----------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            if(debugFlag)
                if(mod(ii,2)==0)
                    clc
                    disp(['iteration: ',num2str(ii),' of ',num2str(nIts)])
                    if(~exist('hp'))
                        hp = plot(20*log10(abs(x)))
                    else
                        set(hp,'YData',20*log10(abs(x)))
                    end
                    title(['res: ',num2str(20*log10(sqrt(rsnew))),' dB'])
                    drawnow;
%                     pause
                end
            end
            
            if useNfft == 1
                Ap = ifft(fftA.*fft(p, nfft));              % Projection of A onto P (How much of A lives on P) -> A*p
            else
                %Ap = ifft(fftA.*fft(p));              % Projection of A onto P (How much of A lives on P) -> A*p
                Ap = ifft(fftA.*fft([p; zeros(nfft-N, 1)]));              % Projection of A onto P (How much of A lives on P) -> A*p
            end
            Ap = Ap(N+[0:M-1]);
            alpha = rsold/(Ap'*Ap);                     % Ratio of Power in Residual / Power of P,  A Matrix Projection (scales P properly so Ax =b)
            x = x + alpha*p;                            % Update Solution (Add this iterations orthogonal guess/basis (P) to existing)
            r = r - alpha*Ap;                           % Project this guess out of residual --> This is negative gradient (towards solution)
            if useNfft == 1
                AhatR = ifft(fftAFC.*fft(r, nfft));         % A'r
            else
                %AhatR = ifft(fftAFC.*fft(r));         % A'r
                AhatR = ifft(fftAFC.*fft([r; zeros(nfft-M, 1)]));         % A'r
            end
            AhatR = AhatR(M + (0:N-1));                 % A'r
            rsnew = sum(abs(AhatR).^2);                 % Computed New Residual (trying to drive this to zero --> means Ax = b)
            
            rOut(ii + 1) = rsnew;
            if sqrt(rsnew) < tol            % Break if pretty darn close
                nIts = ii;
                break;
            end
            
            p = AhatR + rsnew/rsold*p; % Contains all n orthogonal basis vectors (the Ps) (with there powers normalized to power of newest basis vector)
            rsold = rsnew;             % update power number
        end
        
    case 1                                  % This is Block Toeplitz
        if ~isstruct(AStick)                % Even size blocks
            M = size(b, 1);                 % # of rows in A matrix
            nStick = size(AStick, 1);
            N = size(AStick, 1) - M + 1;    % # of delay taps in 1 block of A matrix
            nBlk = size(AStick, 2);         % Number of Blocks
            NTot = N*nBlk;                  % Total Number of taps in A matrix
            if ~exist('x', 'var') || isempty(x)
                x = zeros(N, nBlk);         % Initialize to all zeros
            end
            
            if ~exist('nIts', 'var') || isempty(nIts)
                nIts = NTot;                % Default to Max Rank
            end
            rOut = zeros(1, nIts + 1);

            nfft = 2^(ceil(log2(M + N - 1)));  % Number of Points in FFT
            %nfft = M + N - 1;
            
            
            if useGpu == 1
                AStick = gpuArray(AStick);
                b = gpuArray(b);
                rOut = gpuArray(rOut);
                x = gpuArray(x);
            end
            
            
            AStickFC = conj(flipud(AStick));  % Conj and flip so convolution == xcorr (A')
            if useNfft == 1
                fftA = fft(AStick, nfft);     % (A)
                fftAFC = fft(AStickFC, nfft); % (A')
            else
%                 fftA = fft(AStick);         % (A)
%                 fftAFC = fft(AStickFC);     % (A')
                fftA = fft([AStick; zeros(nfft - nStick, nBlk)]);     % (A)
                fftAFC = fft([AStickFC; zeros(nfft - nStick, nBlk)]); % (A')
            end
            
            
            % % This Computes an XCORR
            %   or = b;
            %   otemp = ifft(fftAFC.*fft(or,nfft)); % xcorr A with b (A'*b)
            %   op = otemp(M+[0:N-1]);              % Grab taps we care about (A'*b)
            %   gamma1 = op'*op;                    % Compute Power b'AA'b
            %
            
            if useNfft == 1
                temp = ifft(fftA.*fft(x, nfft));  % A*x *** For all blocks
            else
                %temp = ifft(fftA.*fft(x));       % A*x *** For all blocks
                temp = ifft(fftA.*fft([x; zeros(nfft-N, nBlk)]));  % A*x *** For all blocks
            end
            r = b - sum(temp(N+[0:M-1], :), 2);  % r = b - A*x , length M *** Sum Across Blocks
            if useNfft == 1
                temp = ifft(fftAFC.*repmat(fft(r,nfft), 1, nBlk));  % xcorr A with r (A'*r)  *** For all Blocks
            else
                %temp = ifft(fftAFC.*repmat(fft(r), 1, nBlk));  % xcorr A with r (A'*r)  *** For all Blocks
                temp = ifft(fftAFC.*repmat(fft([r; zeros(nfft-M, 1)]), 1, nBlk));  % xcorr A with r (A'*r)  *** For all Blocks
            end
            p = temp(M+[0:N-1], :);   % Grab taps we care about (p = A'*r) *** Sum Across Blocks
            rsold =sum(abs(p(:)).^2);
            
            rOut(1) = rsold;  % Hold Residuals
            
            for ii=1:nIts     % Loop over all dimensions
                %-----------------------------------------------------------------------------------------------------------------------------------------------------------------
                %
                % % -------------------- CODE DOES THIS USING FFTs ----------------
                %    %Ap = A*p;                  % Projection of A onto P (How much of A lives on P)
                %    %alpha = rsold/(Ap'*Ap);    % Ratio of Power in Residual / Power of P,  A Matrix Projection (scales P properly so Ax =b)
                %    %x = x + alpha*p;           % Update Solution (Add this iterations orthogonal guess/basis (P) to existing)
                %    %r = r - alpha*Ap;          % Project this guess out of residual --> This is negative gradient (towards solution)
                %    %rsnew = (r'*A)*(A'*r);     % Computed New Residual (trying to drive this to zero --> means Ax = b
                %    %p = A'*r + rsnew/rsold*p;  % Contains all n orthogonal basis vectors (the Ps) (with there powers normalized to power of newest basis vector)
                %    %rsold = rsnew;             % update power number
                %-----------------------------------------------------------------------------------------------------------------------------------------------------------------
                if useNfft == 1
                    Ap = ifft(fftA.*fft(p, nfft));   % Projection of A onto P (How much of A lives on P) -> A*p
                else
                    %Ap = ifft(fftA.*fft(p));        % Projection of A onto P (How much of A lives on P) -> A*p
                    Ap = ifft(fftA.*fft([p; zeros(nfft-N, nBlk)])); % Projection of A onto P (How much of A lives on P) -> A*p
                end
                Ap = sum(Ap(N+[0:M-1], :), 2);
                alpha = rsold/(Ap'*Ap);   % Ratio of Power in Residual / Power of P,  A Matrix Projection (scales P properly so Ax =b)
                x = x + alpha*p;          % Update Solution (Add this iterations orthogonal guess/basis (P) to existing)
                r = r - alpha*Ap;         % Project this guess out of residual --> This is negative gradient (towards solution)
                if useNfft == 1
                    AhatR = ifft(fftAFC.*repmat(fft(r, nfft), 1, nBlk)); % A'r
                else
                    %AhatR = ifft(fftAFC.*repmat(fft(r), 1, nBlk));      % A'r
                    AhatR = ifft(fftAFC.*repmat(fft([r; zeros(nfft-M, 1)]), 1, nBlk));  % A'r
                end
                
                AhatR = AhatR(M + (0:N-1), :); % A'r
                rsnew = sum(abs(AhatR(:)).^2); % Computed New Residual (trying to drive this to zero --> means Ax = b)
                
                rOut(ii + 1) = rsnew;
                if sqrt(rsnew) < tol  % Break if pretty darn close
                    nIts = ii;
                    break;
                end
                
                p = AhatR + rsnew/rsold*p;  % Contains all n orthogonal basis vectors (the Ps) (with there powers normalized to power of newest basis vector)
                rsold = rsnew;              % update power number
            end
        else                           % Uneven Block Sizes (Input as structure)
            
            M = size(b, 1);           % # of rows in A matrix
            nBlk = numel(AStick);     % Number of Blocks
            for ii = 1:nBlk
                AStick(ii).N = size(AStick(ii).AStick, 1) - M + 1; % # of delay taps in 1 block of A matrix
            end
            NTot = sum([AStick(:).N]);  % Total Number of taps in A matrix
            
            if ~exist('x', 'var') || isempty(x)
                for ii = 1:nBlk
                    AStick(ii).x = zeros(AStick(ii).N, 1);
                end
            end
            
            if ~exist('nIts', 'var') || isempty(nIts)
                nIts = NTot;       % Default to Max Rank
            end
            rOut = zeros(1, nIts + 1);
            
            for ii = 1:nBlk
                %AStick(ii).nfft = 2^(ceil(log2(M + AStick(ii).N - 1)));               % Number of Points in FFT
                AStick(ii).nfft = (M + AStick(ii).N - 1);  % Number of Points in FFT
                AStick(ii).AStickFC = conj(flipud(AStick(ii).AStick));  % Conj and flip so convolution == xcorr (A')
                
                if AStick(ii).nfft ~= (M + AStick(ii).N - 1);
                    AStick(ii).fftA = fft(AStick(ii).AStick, AStick(ii).nfft); % (A)
                    AStick(ii).fftAFC = fft(AStick(ii).AStickFC, AStick(ii).nfft); % (A')
                else
                    AStick(ii).fftA = fft(AStick(ii).AStick);     % (A)
                    AStick(ii).fftAFC = fft(AStick(ii).AStickFC); % (A')
                end
            end
            
            %       % This Computes an XCORR
            %         or = b;
            %         otemp = ifft(fftAFC.*fft(or,nfft)); % xcorr A with b (A'*b)
            %         op = otemp(M+[0:N-1]);              % Grab taps we care about (A'*b)
            %         gamma1 = op'*op;                    % Compute Power b'AA'b
            %
            
            r = b;                 % r = b - A*x , length M *** Sum Across Blocks
            for ii = 1:nBlk
                if AStick(ii).nfft ~= (M + AStick(ii).N - 1);
                    temp = ifft(AStick(ii).fftA.*fft(AStick(ii).x, AStick(ii).nfft)); % A*x *** For all blocks
                else
                    temp = ifft(AStick(ii).fftA.*fft(AStick(ii).x)); % A*x *** For all blocks
                end
                r = r - temp(AStick(ii).N+[0:M-1]);  % r = b - A*x , length M *** Sum Across Blocks
            end
            
            for ii = 1:nBlk
                if AStick(ii).nfft ~= (M + AStick(ii).N - 1);
                    temp = ifft(AStick(ii).fftAFC.*fft(r,AStick(ii).nfft)); % xcorr A with r (A'*r)  *** For all Blocks
                else
                    temp = ifft(AStick(ii).fftAFC.*fft(r));    % xcorr A with r (A'*r)  *** For all Blocks
                end
                AStick(ii).p = temp(M+[0:AStick(ii).N-1], :);  % Grab taps we care about (p = A'*r) *** Sum Across Blocks
            end
            rsold = sum(sum(abs([AStick(:).p]).^2));
            
            rOut(1) = rsold;  % Hold Residuals
            
            for ii=1:nIts     % Loop over all dimensions
                %-----------------------------------------------------------------------------------------------------------------------------------------------------------------
                %
                % % -------------------- CODE DOES THIS USING FFTs ----------------
                %    %Ap = A*p;                  % Projection of A onto P (How much of A lives on P)
                %    %alpha = rsold/(Ap'*Ap);    % Ratio of Power in Residual / Power of P,  A Matrix Projection (scales P properly so Ax =b)
                %    %x = x + alpha*p;           % Update Solution (Add this iterations orthogonal guess/basis (P) to existing)
                %    %r = r - alpha*Ap;          % Project this guess out of residual --> This is negative gradient (towards solution)
                %    %rsnew = (r'*A)*(A'*r);     % Computed New Residual (trying to drive this to zero --> means Ax = b
                %    %p = A'*r + rsnew/rsold*p;  % Contains all n orthogonal basis vectors (the Ps) (with there powers normalized to power of newest basis vector)
                %    %rsold = rsnew;             % update power number
                %-----------------------------------------------------------------------------------------------------------------------------------------------------------------
                
                Ap= zeros(M, 1);
                for jj = 1:nBlk
                    if AStick(jj).nfft ~= (M + AStick(jj).N - 1);
                        temp = ifft(AStick(jj).fftA.*fft(AStick(jj).p, AStick(jj).nfft)); % Projection of A onto P (How much of A lives on P) -> A*p
                    else
                        temp = ifft(AStick(jj).fftA.*fft(AStick(jj).p));  % Projection of A onto P (How much of A lives on P) -> A*p
                    end
                    Ap = Ap + temp(AStick(jj).N+[0:M-1]);
                end
                alpha = rsold/(Ap'*Ap);  % Ratio of Power in Residual / Power of P,  A Matrix Projection (scales P properly so Ax =b)
                
                for jj = 1:nBlk
                    AStick(jj).x = AStick(jj).x + alpha*AStick(jj).p; % Update Solution (Add this iterations orthogonal guess/basis (P) to existing)
                end
                r = r - alpha*Ap;  % Project this guess out of residual --> This is negative gradient (towards solution)
                
                for jj = 1:nBlk
                    if AStick(jj).nfft ~= (M + AStick(jj).N - 1)
                        temp = ifft(AStick(jj).fftAFC.*fft(r, AStick(jj).nfft)); % A'r
                    else
                        temp = ifft(AStick(jj).fftAFC.*fft(r)); % A'r
                    end
                    AStick(jj).AhatR = temp(M + (0: AStick(jj).N-1)); % A'r
                end
                rsnew = sum(sum(abs([AStick(:).AhatR]).^2)); % Computed New Residual (trying to drive this to zero --> means Ax = b)
                
                rOut(ii + 1) = rsnew;
                if sqrt(rsnew) < tol % Break if pretty darn close
                    nIts = ii;
                    break;
                end
                for jj = 1:nBlk
                    AStick(jj).p =  AStick(jj).AhatR + rsnew/rsold* AStick(jj).p; % Contains all n orthogonal basis vectors (the Ps) (with there powers normalized to power of newest basis vector)
                end
                rsold = rsnew;  % update power number
            end
            
            for ii = 1:nBlk
                x(ii).x = AStick(ii).x;
            end
        end
end


% % CgOne
% % an implementation of the conjugate gradient algorithm
% %
% % It is used to compute the least squares solution to
% % Ax = b where A is the desired signal and b is the observed signal and x
% % is the solution
% %
% % INPUTS:
% %  dat: observed signal b
% %  ref: desired signal A
% %  alp: iteration limit
% %  tol: residual error limit
% %  M: number of elements in solution
% %  L: number of data points
% %  FFTLen: length of fft used for the matrix vector multiply
% %  iguess: the initial guess
% %
% % OUTPUTS:
% %  x: the answer
% %  k: the number of iterations
% %  BG: norm of the risiduals at each iteration
% %  gammal: norm of data at first iteration
% %
% function [ x,k,BG,gamma1 ] = CgOne( dat,ref,alp,tol,M,L,FFTlen,iguess )
% 
%     refh = conj(fliplr(ref));                 %Desired or A matrix
%     Tfft = fft(ref,FFTlen);                   % FFT of A
%     Thfft = fft(refh,FFTlen);                 % FFT of conj(fliplr(A))
% 
%     r = dat;                                  % b, or observed
%     temp = ifft(Thfft.*fft(r,FFTlen));        % xcorr A with b (A'*b)
%     p = temp(L+[0:M-1]);                      % Grab taps we care about (A'*b)
%     gamma1 = norm(p,2).^2;                    % Compute Power b'AA'b
% 
%     x = iguess;
%     temp = ifft(Tfft.*fft(x, FFTlen));        % b'A?
%     r = dat - temp(M+[0:L-1]);                % b - b'A , length L
% 
%     temp = ifft(Thfft.*fft(r,FFTlen));        % xcorr A with b (A'*b)
%     p = temp(L+[0:M-1]);                      % Grab taps we care about (A'*b)
%     gamma = norm(p,2).^2;
%     BG = zeros(1,alp);
%     BG(1) = gamma;
% 
%     for k=1:alp
% 
%         temp = ifft(Tfft.*fft(p, FFTlen));        % b'A?
%         q = temp(M+[0:L-1]);                      % b'A , length L
%         alpha = gamma./norm(q,2).^2;              % Power of p / power in A
%         x = x + alpha.*p;                         % update weights (length M)
%         r = r - alpha.*q;                         % update residual (length L)
% 
%         temp = ifft(Thfft.*fft(r,FFTlen));        % xcorr A with r (A'*r)
%         s = temp(L+[0:M-1]);                      % Grab taps we care about (A'*r)
%         gammalast = gamma;
%         gamma = norm(s,2).^2;                     % update gamma
%         BG(k) = gamma;
% 
%         if (sqrt(gamma./gamma1) < tol)
%             break;
%         end
% 
%         beta = gamma./gammalast;
%         p = s + beta.*p;
% 
%     end
% 
end