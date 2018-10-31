function result = BM(fun, x) % Balls method
    EPS = 1e-6;
    normA = norm(fun.A);

    grad = fun.Grad(x);
    tetta = acos((sum(-x .* grad)) / norm(x) / norm(grad));
    
    while abs(tetta) > EPS
        c = x - grad / normA;
        a = c' * fun.A * c;
        b = -2 * c' * fun.A * fun.c;
        d = fun.c' * fun.A * fun.c - 1;
        D = b^2 - 4 * a * d;
        t = (-b - sqrt(D)) / 2 / a;
        
        x = t * c;
        grad = fun.Grad(x);
        tetta = acos((sum(-x .* grad)) / norm(x) / norm(grad));
    end
    
    result = x;
end