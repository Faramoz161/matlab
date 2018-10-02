function res = BM(fun, x)
    EPS = 1E-6;

    grad = fun.Grad(x);
    tetta = acos((sum(-x.*grad))/(norm(x)*norm(grad)));
    normA = norm(fun.A);
    
    while(abs(tetta) > EPS)
        c = x - grad/normA;
        a = c.'*fun.A*c;
        b = -2*c.'*fun.A*fun.c;
        d = fun.c.'*fun.A*fun.c - 1;
        D = b*b - 4*a*d;
        t = (-b - sqrt(D))/(2*a);
        
        x = t*c;
        grad = fun.Grad(x);
        tetta = acos((sum(-x.*grad))/(norm(x)*norm(grad)));
    end
    
    res = x;
end