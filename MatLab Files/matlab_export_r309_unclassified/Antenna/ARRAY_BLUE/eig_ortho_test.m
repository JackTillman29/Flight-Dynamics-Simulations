a = rand(3,3);
b = a*a'; % b will have orthonormal eigenvectors

[v,d]=eig(b);

fprintf('b=v*d*inv(v)\n');
disp(v*d*inv(v))
fprintf('b=v*d*v^T\n');
disp(v*d*v')

% this proves that v is orthonormal, e.g. v^-1 = v^T
% should be able to decopmpose the decomposition into a summation series

v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);

d1 = d(1,1);
d2 = d(2,2);
d3 = d(3,3);

sum_result = d1*v1*v1'+d2*v2*v2'+d3*v3*v3';
fprintf('b=sum of lambda*vi*vi''\n');
disp(sum_result);