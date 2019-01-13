function result = PM_1set_Fletcher_Reeves(fun) % Penalty method
    EPS = 1e-4;
    
    n = length(fun.c);
    x = zeros(n, 1);
    r = 1;

    while H(fun, x) > EPS
        x = Fletcher_Reeves(fun, x, r);
        r = r * 2;
    end
    
    result = x;
end

function result = Fletcher_Reeves(fun, initialPoint, r)
    EPS = 1e-5;
    
    x = initialPoint;
    d = - F_grad(fun, x, r);

    n = length(x);    
    iteration = 1;
    while norm(F_grad(fun, x, r)) > EPS
        a = GoldenRatio(fun, x, d, r);
        nextX = x + a * d;
        
        if mod(iteration + 1, n) ~= 0
            b = (norm(F_grad(fun, nextX, r))^2 / norm(F_grad(fun, x, r)))^2;
        else
            b = 0;
        end
        
        d = - F_grad(fun, nextX, r) + b * d;
        iteration = iteration + 1;
        x = nextX;
    end
    
    result = x;
end

function result = GoldenRatio(fun, x_, d_, r)
    EPS = 1e-6;
    LEFT_CONST = (3 - sqrt(5)) / 2;
    RIGHT_CONST = (sqrt(5) - 1) / 2;
    
    a = 0;
    b = 10;
    c = a + LEFT_CONST * (b - a);
    d = a + RIGHT_CONST * (b - a);
    
    fc = F(fun, x_ + c * d_, r);
    fd = F(fun, x_ + d * d_, r);
    
    while b - a > 2 * EPS
        if fc <= fd
            b = d;
            d = c;
            c = a + LEFT_CONST * (b - a);
            
            fd = fc;
            fc = F(fun, x_ + c * d_, r);
        else
            a = c;
            c = d;
            d = a + RIGHT_CONST * (b - a);
            
            fc = fd;
            fd = F(fun, x_ + d * d_, r);
        end
    end
    
    result = (a + b) / 2;
end

function result = F(fun, x, r)
    result = norm(x)^2 + r * H(fun, x);
end

function result = F_grad(fun, x, r)
    result = 2 * (x + r * fun.Val(x) * fun.Grad(x));
end

function result = H(fun, x)
    result = fun.Val(x)^2;
end