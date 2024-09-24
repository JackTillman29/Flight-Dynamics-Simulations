classdef KalmanFilter
    properties
        m_Phi;  % State Transition Matrix
        m_Pk;   % Post Measurement Covariance
        m_Mk;   % Pre Measurement Covariance
        m_Rk;   % Sensor Measurement Covariance
        m_H;    % Mapping Matrix (States -> Measurements)
        m_Qk;   % Process Noise
        m_xhat; % State Estimate
        m_nUpdates = 0; % Counter of number of updates
    end
    methods
        % Constructor
        function obj = KalmanFilter(Phi,H,Rk,Qk)
            obj.m_Phi = Phi;
            obj.m_H = H;
            obj.m_Rk = Rk;
            obj.m_Qk = Qk;
            obj.m_xhat = zeros(size(Phi,1),1);
            obj.m_Pk = 0*(obj.m_xhat * obj.m_xhat');
        end
        
        function obj = InitializeFilter(obj,xhat,Pk)
            if(sum(size(xhat) == size(obj.m_xhat)) == 2)
                obj.m_xhat = xhat;
            else
                error('Wrong size xhat passed in');
            end
            if(sum(size(Pk) == size(obj.m_Pk)) == 2)
                obj.m_Pk   = Pk;
            else
                error('Wrong size Pk passed in');
            end
        end
        
        function obj = UpdateMeasurement(obj,xm)
            
        end
        
        % Static Functions
        % Compute the discrete state transition matrix
        function Phi = ComputePhi(filter,F,dt)
            Phi = expm(F.*dt);
        end
        
        % Compute the discrete process noise
        function [Qk] = ComputeNumericQk(filter,F,Q,Ts,nPts)
            
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
        
    end
end