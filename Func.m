classdef Func    
    properties(SetAccess = public)
        A;
        c;
    end
    
    methods
        function obj = Func(A, c)
            obj.A = A;
            obj.c = c;
        end
        
        function res = Val(this, x)
           res = (x - this.c).'*this.A*(x - this.c) - 1;
        end
        function res = Grad(this, x)
            res = 2*this.A*(x - this.c);
        end
        function res = H(this)
            res = 2*this.A;
        end 
    end
end

