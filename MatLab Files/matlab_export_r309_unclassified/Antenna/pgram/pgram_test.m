% all units are in meters
% array lies (initially) in x-y plane

[X,Y] = meshgrid(0:0.01:1,0:0.01:1);

close all;
figure;
plot(X(:),Y(:),'.');

lineA = line([0 0.3],[0 1]);
lineB = line([0 0.7],[0 0]);
lineC = line([1 0.7],[1 0]);
lineD = line([1 0.3],[1 1]);

set([lineA lineB lineC lineD],'Color','r');

% compute cross product vectors
u = [X(:) Y(:) 0*X(:)];
v = [-X(:)+1 -Y(:)+1 0*X(:)];
A = [0.3 1 0];
B = [0.7 0 0];
C = [-0.3 -1 0];
D = [-0.7 0 0];


uc = zeros(length(u),4);

for k = 1 : length(u)
    cpA = cross(u(k,:), A);
    cpB = cross(u(k,:), B);
    cpC = cross(v(k,:), C);
    cpD = cross(v(k,:), D);
    uc(k,:) = [ ...
        cpA(3) cpB(3) cpC(3) cpD(3) ];
    
end

good = 0 * X;
p = find( ...
    uc(:,1) > 0 & ...
    uc(:,2) < 0 & ...
    uc(:,3) < 0 & ...
    uc(:,4) > 0);
good(p) = 1;

figure;
scatter(X(:),Y(:),108,good(:),'filled','sq');
colorbar;

title('test 4');