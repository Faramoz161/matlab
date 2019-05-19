function x = PM_projection(fun, initialPoint)
    % F = norm(x)^2 + r * fun.Val(x)^2
    
    EPS = 1e-6;
    
    x = initialPoint;
    r = 1000;

    while norm(Psi(fun, x)) > EPS
        
        while norm(F_grad(fun, x, r)) > EPS / 10
            x = x - F_H(fun, x, r) \ F_grad(fun, x, r);
        end

        grad = fun.Grad(x);
        x = x - grad * fun.Val(x) / norm(grad)^2;
        r = r * 10;
    end
end

function result = F_grad(fun, x, r)
    result = 2 * (x + r * fun.Val(x) * fun.Grad(x));
end

function result = F_H(fun, x, r)
    grad = fun.Grad(x);
    result = 2 * (eye(length(x)) + r * (grad * grad' + fun.Val(x) * fun.H()));
end