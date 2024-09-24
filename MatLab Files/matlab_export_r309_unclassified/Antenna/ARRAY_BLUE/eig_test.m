close all;
clc;

% Create fake eigendecomp
temp = rand(4,4);
A = temp*temp'; % orthonormal

[V,D]=eig(A);
nz = 0.1;

J1 = 10;
J2 = 5;

D(1,1) = J1+nz;
D(2,2) = J2+nz;
D(3,3) =    nz;
D(4,4) =    nz;

% Create new A
A = V*D*V';

%% =========================================================================
M = A;
%[V,D] = eig(M);  % Decompose

theAnswer = nz * inv(M)

% Form corrlation matrix using noise outside the summing junction with
% known jammer eigenvalues
M1 = nz*eye(4) + ...
    J1 * V(:,1)*V(:,1)' + ...
    J2 * V(:,2)*V(:,2)';

% Form correlation matrix using noise outside summing junction with noise
% subtracted from computed eigenvalues: EV = Jam - Noise
M2 = nz*eye(4) + ...
    (D(1,1)-nz) * V(:,1)*V(:,1)' + ...
    (D(2,2)-nz) * V(:,2)*V(:,2)';

% Same. Now invert M2 x noise
M2inv = inv(M2) .* nz

% Now do my derivation
M2inv_kds = ...
    1./(1+(D(1,1)-nz)./nz) * V(:,1)*V(:,1)' + ...
    1./(1+(D(2,2)-nz)./nz) * V(:,2)*V(:,2)' + ...
    1./(1+(D(3,3)-nz)./nz) * V(:,3)*V(:,3)'+ ...
    1./(1+(D(4,4)-nz)./nz) * V(:,4)*V(:,4)'

% Now do Farina derivation
M2inv_farina = ...
    (1 - (D(1,1)-nz)./D(1,1)) * V(:,1)*V(:,1)' + ...
    (1 - (D(2,2)-nz)./D(2,2)) * V(:,2)*V(:,2)' + ...
    (1 - (D(3,3)-nz)./D(3,3)) * V(:,3)*V(:,3)' + ...
    (1 - (D(4,4)-nz)./D(4,4)) * V(:,4)*V(:,4)'

% Spin on kds 
N = 4;
sumMatrix = 0;
for k = 1 : N
    sumMatrix = sumMatrix + ...
        1./(1+(D(k,k)-nz)./nz) * V(:,k)*V(:,k)';
end
disp('original');
sumMatrix

sumMatrix = 0;
for k = 1 : N
    sumMatrix = sumMatrix + ...
        (nz./D(k,k)) * V(:,k)*V(:,k)';
end
disp('factored');
sumMatrix


% Create a complex array
R = [ ...
    1.0+1i*5 0.2+1i*3 0.5-1i*9
    0.2-1i*3 2.0+1i*4 0.4-1i*3
    0.5+1i*9 0.4+1i*3 3.0+1i*2];