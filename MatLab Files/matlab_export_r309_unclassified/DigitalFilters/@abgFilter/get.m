function value = get(r,field)
if ( nargin == 2 )
    value = eval(['r.' field]);
else
    disp('valid options are: ');
    nf = fields(r);
    for k = 1:length(nf)
        disp(nf{k});
    end
    disp(' ');
end
end