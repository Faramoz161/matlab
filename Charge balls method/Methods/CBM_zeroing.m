function result = CBM_zeroing(fun, initialPoint) % Charged balls method with zeroing speed
    EPS = 1e-6;
    EPS1 = 1e-8;
    EPS2 = 1e-2;
    DELTA = 0.1;
    p1 = 100;
    p2 = 2;
    
    x = initialPoint;
    z = DELTA * (p1 * Psi(fun, x) - Hi(fun, x, zeros(length(x), 1)));
    
    while norm(Psi(fun, x)) > EPS
        xPrev = x;
        zPrev = z;

        while 1
            xPrevMod = xPrev + DELTA * zPrev;
            grad = fun.Grad(xPrevMod);
            x = xPrevMod - grad * fun.Val(xPrevMod) / (norm(grad)^2);
            
            condition = norm(x - xPrevMod) / norm(xPrevMod);
            if condition < EPS1
                 DELTA = DELTA * 2;
            elseif condition > EPS2
                DELTA = DELTA / 2;
            else
                break;
            end
        end
        
        z = zPrev + DELTA * (p1 * Psi(fun, xPrev) - p2 * zPrev - Hi(fun, xPrev, zPrev));
        
        if norm(x) - norm(xPrev) > 0.1
            z = DELTA * (p1 * Psi(fun, x) - Hi(fun, x, zeros(length(x), 1)));
        end
    end
    
    result = x;
end

function result = Psi(fun, x)
    gr = fun.Grad(x);
    result = (gr * x' * gr / norm(gr)^2 - x) / norm(x)^3;
end

function result = Hi(fun, x, z)
    gr = fun.Grad(x);
    result = gr * z' * fun.H() * z / norm(gr)^2;
end