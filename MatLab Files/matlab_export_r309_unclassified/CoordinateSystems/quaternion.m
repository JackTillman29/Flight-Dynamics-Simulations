classdef quaternion
    properties
        q = [ 0 0 0 0 ];
    end
    methods
        function obj = quaternion(q)
            obj.q = q;
        end
        function display(this)
            disp([this.q(1) this.q(2:4)]');
        end
        function boolean = eq(this,q)
            if(strcmp(class(q),'quaternion'))
                boolean = this.q(1) == q.q(1) && ...
                    this.q(2) == q.q(2) && ...
                    this.q(3) == q.q(3) && ...
                    this.q(4) == q.q(4);
            else
                error(['quaternion == operator undefined for arguments of type ' class(q)]);
            end
        end
        % add one quaternion to another
        function obj = plus(this,q)
            if(strcmp(class(q),'quaternion'))
                obj = quaternion(this.q+q.q);
            else
                error(['quaternion + operator undefined for arguments of type ' class(q)]);
            end
        end

        % subtract one quaternion from another
        function obj = minus(this,q)
            if(strcmp(class(q),'quaternion'))
                obj = quaternion(this.q-q.q);
            else
                error(['quaternion + operator undefined for arguments of type ' class(q)]);
            end
        end

        % negate a quaternion
        function obj = uminus(this)
            obj = quaternion(-this.q);
        end

        % scalar times a quaternion
        function obj = times(a,b)
            if(strcmp(class(a),'quaternion'))
                quat = a;
                scalar = b;
            else
                quat = b;
                scalar = a;
            end
            obj = quaternion(scalar * quat.q);
        end
        
        % quaternion times a quaternion
        function obj = mtimes(this,q)
            if(strcmp(class(q),'quaternion'))
                P = [ ...
                this.q(1)      -this.q(2)   -this.q(3) -this.q(4);
                this.q(2)  this.q(1)        -this.q(4)  this.q(3);
                this.q(3)  this.q(4)    this.q(1)      -this.q(2);
                this.q(4) -this.q(3)    this.q(2)  this.q(1) ];
            
                r = P * [q.q(1) q.q(2) q.q(3) q.q(4)]';
                
                obj = quaternion(r');
            else
                error(['quaternion * operator undefined for arguments of type ' class(q)]);
            end
        end
        
        % quaternion conjugate
        function obj = conj(this)
            obj = quaternion([this.q(1) -this.q(2:4)]);
        end
        
        % norm
        function qnorm = norm(this)
            v = conj(this) * this;
            qnorm = sqrt(v.q(1));
        end
        
        % inverse
        function obj = inv(this)
            obj = conj(this) .* (1/norm(this)^2);
        end
        
        % convert to rotation matrix R^3
        function R = toRotationMatrix(this)
             R = [ ...
                 2*this.q(1)^2-1+2*this.q(2)^2 2*this.q(2)*this.q(3)-2*this.q(1)*this.q(4) 2*this.q(2)*this.q(4)+2*this.q(1)*this.q(3);
                 2*this.q(2)*this.q(3)+2*this.q(1)*this.q(4) 2*this.q(1)^2-1+2*this.q(3)^2 2*this.q(3)*this.q(4)-2*this.q(1)*this.q(2);
                 2*this.q(2)*this.q(4)-2*this.q(1)*this.q(3) 2*this.q(3)*this.q(4)+2*this.q(1)*this.q(2) 2*this.q(1)^2-1+2*this.q(4)^2];
             R = R ./ norm(R);
        end
        
        function 
        
    end
end