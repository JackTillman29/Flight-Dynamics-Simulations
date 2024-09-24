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
        xmeas;
        Kk;
        MeasEqn;        % returns measurements as function of states
    end
    
    methods
        function obj = ExtKalmanFilterH(nStates,nMeasurements,hFcnHandle,mFcnHandle)
            obj.xhat = zeros(nStates,1);
            obj.xmeas= zeros(nMeasurements,1);
            obj.Kk   = zeros(nStates,nMeasurements);
            obj.I    = zeros(nStates,nStates);
            obj.Pk   = zeros(nStates,nStates);
            obj.Qk   = zeros(nStates,nStates);
            obj.Phi  = zeros(nStates,nStates);
            obj.Mk   = zeros(nStates,nStates);
            obj.Gk   = zeros(nStates,nStates);
            obj.H    = hFcnHandle;
            obj.MeasEqn = mFcnHandle;
        end
        
        function this = UpdateMeasurement(this,xmeas)
            this.xmeas = xmeas;
            this.RunFilter();
        end
        
        function this = RunFilter(this)
            H  = this.H(this.xhat);           % Compute H matrix this timestep
            this.Mk = this.Phi * this.Pk * this.Phi' + this.Qk;
            this.Kk = this.Mk * H' * inv(H * this.Mk * H' + this.Rk);
            this.Pk = (this.I - this.Kk * H) * this.Mk;
            xbar = this.Phi * this.xhat + this.Gk;
            residual = this.xmeas - this.MeasEqn(xbar);
            this.xhat = xbar + this.Kk * residual;
        end
        
        
    end
    
end

