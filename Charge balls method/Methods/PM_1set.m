function result = PM_1set(fun) 
    % Penalty method
    
    EPS = 1e-6;
    
    n = length(fun.c);
    x = zeros(n, 1);
    r = 1;

    while H(fun, x) > EPS
        
        while norm(F_grad(fun, x, r)) > EPS / 10
            x = x - F_H(fun, x, r) \ F_grad(fun, x, r);
        end

        r = r * 10;
    end
    
    result = x;
end

function result = F_grad(fun, x, r)
    result = 2 * (x + r * fun.Val(x) * fun.Grad(x));
end

function result = F_H(fun, x, r)
    grad = fun.Grad(x);
    result = 2 * (eye(length(x)) + r * (grad * grad' + fun.Val(x) * fun.H()));
end

function result = H(fun, x)
    result = fun.Val(x)^2;
end