classdef Func    
    properties(SetAccess = private)
        A;
        c;
    end
    
    methods
        function obj = Func(A, c)
            obj.A = A;
            obj.c = c;
        end
        
        function result = Val(this, x)
           result = (x - this.c)' * this.A * (x - this.c) - 1;
        end
        
        function result = Grad(this, x)
            result = 2 * this.A * (x - this.c);
        end
        
        function result = H(this)
            result = 2 * this.A;
        end
        
        function result = Psi(this, x)
            grad = this.Grad(x);
            result = (grad * x' * grad / norm(grad)^2 - x) / norm(x)^3;
        end
        
        function result = Hi(this, x, z)
            grad = this.Grad(x);
            result = grad * z' * this.H() * z / norm(grad)^2;
        end
    end
end

