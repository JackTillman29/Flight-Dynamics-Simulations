function [Qk] = ComputeNumericQk(F,Q,Ts,nPts)

    dt = Ts / (nPts-1);
    temp = zeros(nPts,size(F,1)*size(F,2));
    t = 0:dt:Ts;

    for k = 1:nPts
        Phi_k = expm(F*t(k));
        T     = Phi_k * Q * Phi_k';
        temp(k,:) = T(1:end);
    end

    yk = temp(2:end,:);
    ykm1 = temp(1:(end-1),:);

    % trapezoidal integration
    A = cumsum([zeros(size(temp,2)); ykm1+yk]) * dt / 2;

    Qk = reshape(A(end,:) - A(1,:) ,[size(F,1) size(F,2)] );

end