function MixProducts2_Plotter(A,F,P) 
% each frequency will generate two products when mixed with another
% frequency. Note: P is phase ref, and is either 0 or 1 (out of phase)

hLinePos = [];
hLineNeg = [];

for k = 1:length(A)
    if(P(k) == 1)
        hLinePos = [hLinePos line([F(k) F(k)],[0 -A(k)])];
    else
        hLineNeg = [hLineNeg line([F(k) F(k)],[0  A(k)])];
    end
end

set(hLinePos,'LineWidth',2,'Color','b');
set(hLineNeg,'LineWidth',2,'Color','r');

end