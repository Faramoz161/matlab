classdef Func    
    
    properties(SetAccess = public)
        A;
        c;
    end
    
    methods
        
        % Constructor
        function obj = Func(A, c)
            obj.A = A;
            obj.c = c;
        end

        % Calculation f(x)
        function v = Val(this, x)
           v = (x - this.c)' * this.A * (x - this.c) - 1;
        end
        
        % Gradient calculation
        function gr = Grad(this, x)
            gr = 2 * this.A * (x - this.c);
        end
        
        % Calculation of the matrix of second derivatives
        function h = H(this)
            h = 2 * this.A;
        end
        
    end
end

