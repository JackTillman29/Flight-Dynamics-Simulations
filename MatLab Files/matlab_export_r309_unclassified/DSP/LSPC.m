function x = LSPC(b,kernel,AStick)
%  MAY NEED TO CIRCULARLY SHIFT OR REVERSE THE OUTPUT (???)
%  INPUTS: 
%     b      = received signal to process with LSPC
%     kernel = [Nx1 or 1xN]
%  OUTPUT:
%     x      = range profile (point targets are essentially impulse
%                             functions convolved with the kernel/pulse so
%                             LSPC attempts
%
%  Least-Squares Pulse Compression solver
%     Solves:
%       A x = b
%       x = inv(A^H * A) * A^H * b
%
%  For large b and/or kernel, forming the A matrix and solving 
%  inv(A^H * A) * A^H can be prohibitively expensive (in terms of memory 
%  and ops), so it is solved iteratively using
%     conjugate gradient least squares method (CG)
%       *** depends on cgToep (conjugate gradient method to solve:
%
%  Application to radar received signals with a known transmitted pulse:
%    Point targets are essentially impulse functions, and the received 
%    signal is the convolution of these impulse functions with the 
%    kernel/pulse. 
%    In matched filtering pulse compression, the goal is to maximize SNR at
%       the expense of introducing range sidelobes.
%    In LSPC, a form of mismatched filtering (not optimal in SNR sense),
%       the goal is to recover the point targets' impulse functions.


% put inputs into column vectors
if(size(b,2)~=1)
    y_dim = 'row';
    b = b.';
else
    y_dim = 'col';
end
if(size(kernel,2)~=1)
    kernel = kernel.';
end
    
N = length(kernel);
M = length(b);
L = M - N + 1;

%==========================================================================
% INTERFACING TO cgToep
%  x = inv(A^H * A) * A^H * b
%       AStick == A^H ???    (H == hermitian = conj(A'))
%       A^H_(M+N-1xM) * b_Mx1
%
%   M = length(b) = length(y)
%   N = length(kernel)
%       kernel goes into "AStick"
%            hence M+N-1x1 = length(y)+length(kernel)-1 x 1
%
% looking at docs from study, AStick "outer edge" is
%    -<<------------------1st element
%   |          size M
%   |    [ a1 0  0  ...   0  ]    a = [a1 a2 ... aN] is length-N kernel
%   |    [ a2 a1 0  ...   0  ]
%   |    [ a3 a2 a1 ...   0  ]
%   |    [ .  .  .        .  ]  size N+(L-1) = N+(
%   |    [ .  .  .        a0 ]
%   |    [ .  .  .        .  ]
%   V    [ 0  0  0  ...   aN ]
%  (last element)
%   AStick = [zeros(M-1) a1 a2 a3 ... aN zeros(M-1)]

if nargin <3
    AStick = [zeros(M-1,1); conj(kernel); zeros(L-N,1)];
end
% AStick = [zeros(L-N,1); conj(kernel); zeros(M-1,1)];

% BAD
% AStick = [zeros(L-N,1); conj(kernel); zeros(M-1,1)];
% GOOD? BUT NOT ALIGNED TO MATCHED FILTER...
% AStick = [zeros(M-1,1); kernel; zeros(L-N,1)];
% BAD
% AStick = [zeros(L-N,1); kernel; zeros(M-1,1)];

useNFFT  = 0;
toepType = 0; % 0: reg Toeplitz, 1: block Toeplitz
nIterations = 300;
[x,nIts,rOut] = cgToep(AStick, b(1:end), [], ...
                       toepType,[],nIterations,useNFFT);

% flip x to align with true amb range
%    x is ALWAYS a column vector at this point.
x = flipud(x);

% output of cgToep is column vector so form x into shape of y
if(strcmpi(y_dim,'row'))
    x = x.';
end

% % DIRECT METHOD (form Toeplitz S matrix and solve):
% % disp(['N = ',num2str(N)])
% % disp(['L = ',num2str(L)])
% % disp(['S will be Nx(N+L-1) = ',num2str(N),'x',num2str(N+L-1)])
% %  form S using sparse matrix for speedup
% S = sparse( toeplitz([kernel; zeros(L-1,1)],[kernel(1) zeros(1,L-1)]) );
% % S = sparse( toeplitz([kernel; zeros(L-1,1)],[kernel(1) zeros(1,M-1)]) );
% x = inv(S'*S) * (S' * y);
% % x = pinv(S)*y; % will not work with sparse matrix because of call to svd, 
% %                % which does not support sparse matrices (does say to use
% %                % svds for sparse matrix svd, but that would mean modifying 
% %                % pinv)




end