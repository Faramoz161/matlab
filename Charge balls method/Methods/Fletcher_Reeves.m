function result = Fletcher_Reeves(fun, initialPoint)
    EPS = 1e-6;
    
    x = initialPoint;
    d = -F_grad(fun, x);
    
    n = length(x);    
    iteration = 1;
    while norm(Psi(fun, x)) > EPS
        a = GoldenRatio(fun, x, d);
        nextX = x + a * d;
        
        if mod(iteration + 1, n) ~= 0
            b = (norm(F_grad(fun, nextX)) / norm(F_grad(fun, x)))^2;
        else
            b = 0;
        end
        
        d = -F_grad(fun, nextX) + b * d;
        iteration = iteration + 1;
        
        grad = fun.Grad(nextX);
        x = nextX - grad * fun.Val(nextX) / norm(grad)^2;
    end
    
    result = x;
end

function result = GoldenRatio(fun, x_, d_)
    EPS = 1e-6;
    LEFT_CONST = (3 - sqrt(5)) / 2;
    RIGHT_CONST = (sqrt(5) - 1) / 2;
    
    a = 0;
    b = 1;
    c = a + LEFT_CONST * (b - a);
    d = a + RIGHT_CONST * (b - a);
    
    fc = F(fun, x_ + c * d_);
    fd = F(fun, x_ + d * d_);
    
    while b - a > 2 * EPS
        if fc <= fd
            b = d;
            d = c;
            c = a + LEFT_CONST * (b - a);
            
            fd = fc;
            fc = F(fun, x_ + c * d_);
        else
            a = c;
            c = d;
            d = a + RIGHT_CONST * (b - a);
            
            fc = fd;
            fd = F(fun, x_ + d * d_);
        end
    end
    
    result = (a * b) / 2;
end

function result = Psi(fun, x)
    gr = fun.Grad(x);
    result = (gr * x' * gr / norm(gr)^2 - x) / norm(x)^3;
end

function result = F(fun, x)
    result = norm(x)^2 + fun.Val(x)^2;
end

function result = F_grad(fun, x)
    result = 2 * (x + fun.Val(x) * fun.Grad(x));
end

function result = F_H(fun, x)
    grad = fun.Grad(x);
    result = 2 * (eye(length(x)) + grad * grad' + fun.Val(x) * fun.H());
end