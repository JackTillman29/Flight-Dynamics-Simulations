function [Amix,Fmix] = MixProducts(A1,F1,A2,F2) 
% each frequency will generate two products when mixed with another
% frequency
nProducts = 2*length(F1) * length(F2);

if(isempty(A1))
    A1 = 0*F1 + 1;
end

if(isempty(A2))
    A2 = 0*F2 + 1;
end

Amix = zeros(1,nProducts);
Fmix = zeros(1,nProducts);

iProduct = 1;
for k1 = 1:length(F1)
    for k2 = 1:length(F2)
        Fmix(iProduct) = abs(F1(k1)+F2(k2));
        Amix(iProduct) = (A1(k1)*A2(k2))/2;
        iProduct = iProduct + 1;
        Fmix(iProduct) = abs(F1(k1)-F2(k2));
        Amix(iProduct) = (A1(k1)*A2(k2))/2;
        iProduct = iProduct + 1;
    end
end