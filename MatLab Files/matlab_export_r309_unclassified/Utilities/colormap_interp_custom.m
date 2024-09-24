function cmnew = colormap_interp_custom(cm,N)
% colormap_interp_custom(cm,N)
%    cm = Nx3 matrix of RGB values
%    N  = number of interpolation points between RGB values in cm

f1 = @(N,x,y) linspace(x,y,N).';
cmnew = [];
for i = 1:(size(cm,1)-1)
    cm1 = [f1(N,cm(i,1),cm(i+1,1)), ...
           f1(N,cm(i,2),cm(i+1,2)), ...
           f1(N,cm(i,3),cm(i+1,3))];
    cmnew = [cmnew; cm1];
end

end