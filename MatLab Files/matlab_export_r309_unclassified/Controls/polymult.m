function polyout = polymult(poly1,poly2,trim)
if(exist('trim')~=1)
    trim = 0;
elseif(strcmp(trim,'trim'))
    trim = 1;
else
    trim = 0;
end
% zero pad inputs so everyone is the same size
iMaxLength = max([length(poly1) length(poly2)]);

% size deficiency arrays
sda = [ ...
    iMaxLength - length(poly1) ...
    iMaxLength - length(poly2)];

if(sda(1) ~= 0)
    poly1 = [zeros(1,sda(1)) poly1];
end
if(sda(2) ~= 0)
    poly2 = [zeros(1,sda(2)) poly2];
end

% sizes are all okay, now apply the rules


poly_out = zeros(1,2*iMaxLength-1);

% foil numerator
% s^3 s^2 s^1 s^0      s^3 s^2 s^1 s^0
% 
OL = (iMaxLength:-1:1);
for c1 = 1:iMaxLength
    for c2 = 1:iMaxLength
        cout = OL(c1)+OL(c2)-1;
        poly_out(cout) = poly_out(cout) + poly1(c1)*poly2(c2);
    end
end
   
polyout = poly_out(end:-1:1);
if(trim)
nidx = find(polyout ~= 0,1,'first');
polyout = polyout(nidx:end);
end