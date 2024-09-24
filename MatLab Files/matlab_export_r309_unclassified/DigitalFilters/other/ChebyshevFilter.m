function [num,den] = ChebyshevFilter(FC,LH,PR,NP)
% This code is based on "Scientists & Engineers Guide to DSP by Steven
% Smith". Page 340-341.
% function [A,B] = ChebyshevFilter(FC,LH,PR,NP)
% A are the numerator tapped delay gains
% B are the denominator tapped delay gains
% FC is the fraction of the discrete sample frequency to break (i.e. 0.2 on
% a 1kHz system would be 200Hz.
% LH is the low pass (0) or high pass (1) flag. 
% PR is percent ripple. Limited between 0 and 29 (0% equates to
% butterworth)
% NP is the number of poles (even value)
% 
% Coefficient Convention: A0 is applied to Input(k)
%                         A1 is applied to Input(k-1)
%                         A2 is applied to Input(k-2)
%                         B1 is applied to Output(k-1)
%                         B2 is applied to Output(k-2)
% 
%  Y(k)   A0*z^0 + A1*z^-1 + A2*z^-2
%  ---- = --------------------------
%  U(k)   1      - B1*z^-1 - B2*z^-2

if(NP > 20)
    error('Max poles is 20. Do you really need more than that?');
end


% cascade arrays
A  = zeros(1,23);
B  = zeros(1,23);
TA = zeros(1,23);
TB = zeros(1,23);

A(2+1) = 1;
B(2+1) = 1;


% loop over each pole pair
for P = 1 : (NP/2)
    %% GOSUB 1000 Subroutine
    
    % compute pole location on the unit circle
    RP = -cos(pi/(NP*2)+(P-1)*pi/NP);
    IP = +sin(pi/(NP*2)+(P-1)*pi/NP);
    
    
    if ( PR ~= 0.0 ) % if passband ripple
        % warp from circle to ellipse
        ES = sqrt( ( 100 / (100-PR))^2 - 1 );
        VX = (1/NP) * log( ( 1/ES) + sqrt( ( 1/ES^2) + 1 ) );
        KX = (1/NP) * log( ( 1/ES) + sqrt( ( 1/ES^2) - 1 ) );
        KX = (exp(KX) + exp(-KX)) / 2;
        RP = RP * ( ( exp(VX) - exp(-VX) ) / 2 ) / KX;
        IP = IP * ( ( exp(VX) + exp(-VX) ) / 2 ) / KX;
    end
    
    T = 2 * tan(1/2);
    W = 2*pi*FC;
    M = RP^2 + IP^2;
    D = 4.0 - 4.0*RP*T + M*T^2;
    X0 = T^2/D;
    X1 = 2*T^2/D;
    X2 = T^2/D;
    Y1 = (8-2*M*T^2)/D;
    Y2 = (-4 -4*RP*T - M*T^2)/D;
    
    % low pass to low pass or highpass
    if(LH == 1)
        K = -cos( W/2 + 1/2) / cos(W/2 - 1/2);
    else
        K = +sin(-W/2 + 1/2) / sin(W/2 + 1/2);
    end
    
    D = 1 + Y1*K - Y2*K^2;
    A0 = (X0 - X1*K + X2*K^2)/D;
    A1 = (-2*X0*K + X1 + X1*K^2 - 2*X2*K)/D;
    A2 = (X0*K^2 -X1*K + X2)/D;
    B1 = (2*K + Y1 + Y1*K^2 - 2*Y2*K)/D;
    B2 = (-K^2 -Y1*K+Y2)/D;
    
    if(LH==1)
        A1 = -A1;
        B1 = -B1;
    end 
    % end of SUB1000
    %% continue
    for i = 0:22
        TA(i+1) = A(i+1);
        TB(i+1) = B(i+1);
    end
    
    for i = 2:22
        A(i+1) = A0*TA(i+1) + A1*TA(i) + A2*TA(i-1);
        B(i+1) =    TB(i+1) - B1*TB(i) - B2*TB(i-1);
    end
    
end % next pole pair

B(2+1) = 0;
for i = 0 : 20
    A(i+1) =  A(i+2+1);
    B(i+1) = -B(i+2+1);
end

%% normalize the gain
SA = 0;
SB = 0;

for i = 0:20
    if(LH==1)
        SA = SA + A(i+1) * ((-1)^(i));
        SB = SB + B(i+1) * ((-1)^(i));
    else
        SA = SA + A(i+1);
        SB = SB + B(i+1);
    end
end

GAIN = SA / ( 1 - SB );

for i = 0 : 20
    A(i+1) = A(i+1) / GAIN;
end

% reduce data sizes for output
A = A(1:(NP+1));
B = B(1:(NP+1));
num = A;
den = -B;
den(1) = 1;

pole_mag = abs(roots(den));
if( sum( pole_mag > 1 ) )
    warning('Filter is unstable!');
end

end