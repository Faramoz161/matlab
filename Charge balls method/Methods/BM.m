function result = BM(fun, initialPoint) 
    % Balls method

    EPS = 1e-6;
    normA = norm(fun.A);
    
    x = initialPoint;
    grad = fun.Grad(x);
    tetta = acos((-x)' * grad / norm(x) / norm(grad));
    
    d = fun.c' * fun.A * fun.c - 1;
    
    while abs(tetta) > EPS
        c = x - grad / normA;
        cA = c' * fun.A;
        
        a = cA * c;
        b = -2 * cA * fun.c;
        
        D = b^2 - 4 * a * d;
        t = (-b - sqrt(D)) / (2 * a);
        
        x = t * c;
        grad = fun.Grad(x);
        tetta = acos((-x)' * grad / norm(x) / norm(grad));
    end
    
    result = x;
end