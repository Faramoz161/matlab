function x = CBM(fun, initialPoint)
    EPS = 1e-6;
    EPS1 = 1e-8;
    EPS2 = 1e-2;
    DELTA = 1;
    p1 = 20;
    p2 = 1;
    
    x = initialPoint;
    z = DELTA * (p1 * fun.Psi(x) - fun.Hi(x, zeros(length(x), 1)));
    
    while norm(fun.Psi(x)) > EPS
        xPrev = x;
        zPrev = z;
        
        while 1
            xPrevMod = xPrev + DELTA * zPrev;
            grad = fun.Grad(xPrevMod);
            x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
            
            condition = norm(x - xPrevMod) / norm(xPrevMod);
            if condition < EPS1
                DELTA = DELTA * 2;
            elseif condition > EPS2
                DELTA = DELTA / 2;
            else
                break;
            end
        end
        
        z = zPrev + DELTA * (p1 * fun.Psi(xPrev) - p2 * zPrev - fun.Hi(xPrev, zPrev));
    end
end