classdef ExtKalmanFilterH
    properties
        xhat;
        I;
        Pk;
        Rk;
        Qk;
        Phi;
        Mk;
        H;
        Gk;
    end
    
    methods
        function obj = ExtKalmanFilterH(nStates,hFcnHandle)
            obj.xhat = zeros(nStates,1);
            obj.I    = zeros(nStates,nStates);
            obj.Pk   = zeros(nStates,nStates);
            obj.Qk   = zeros(nStates,nStates);
            obj.Phi  = zeros(nStates,nStates);
            obj.Mk   = zeros(nStates,nStates);
            obj.Gk   = zeros(nStates,nStates);
            obj.H    = hFcnHandle;
        end
    end
    
end

