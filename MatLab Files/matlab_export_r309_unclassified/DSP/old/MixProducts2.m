function [Amix,Fmix,Pmix] = MixProducts2(A1,F1,P1,A2,F2,P2) 
% each frequency will generate two products when mixed with another
% frequency. Note: P is phase ref, and is either 0 or 1 (out of phase)
nProducts = 2*length(F1) * length(F2);

if(isempty(A1))
    A1 = 0*F1 + 1;
end

if(isempty(A2))
    A2 = 0*F2 + 1;
end

if(isempty(P1))
    P1 = 0*F1 + 1;  % zero, 0 deg
end

if(isempty(P2))
    P2 = 0*F2 + 1;  % zero, 0 deg
end

Amix = zeros(1,nProducts);
Fmix = zeros(1,nProducts);
Pmix = zeros(1,nProducts);

iProduct = 1;
for k1 = 1:length(F1)
    for k2 = 1:length(F2)
        Fmix(iProduct) = F1(k1)+F2(k2);
        Amix(iProduct) = (A1(k1)*A2(k2))/2;
        Pmix(iProduct) = mod(P1(k1)+P2(k2),2);
        if ( Fmix(iProduct) < 0.0 )
            Pmix(iProduct) = 1-Pmix(iProduct); % phase sense change
            Fmix(iProduct) = -1.0 * Fmix(iProduct);
        end
        iProduct = iProduct + 1;
        
        Fmix(iProduct) = F1(k1)-F2(k2);
        Amix(iProduct) = (A1(k1)*A2(k2))/2;
        Pmix(iProduct) = mod(P1(k1)-P2(k2),2);
        if ( Fmix(iProduct) < 0.0 )
            Pmix(iProduct) = 1-Pmix(iProduct); % phase sense change
            Fmix(iProduct) = -1.0 * Fmix(iProduct);
        end
        iProduct = iProduct + 1;
    end
end