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
        function result = Val(this, x)
           result = (x - this.c)' * this.A * (x - this.c) - 1;
        end
        
        % Gradient calculation
        function result = Grad(this, x)
            result = 2 * this.A * (x - this.c);
        end
        
        % Calculation of the matrix of second derivatives
        function result = H(this)
            result = 2 * this.A;
        end
        
    end
end

