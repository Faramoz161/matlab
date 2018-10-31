function result = CBM_constStep(fun, initialPoint) % Charged balls method with constant step
    EPS = 1e-6;
    DELTA = 0.3;
    p1 = 10;
    p2 = 5;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    
    while norm(Psi(fun, x)) > EPS
        xPrev = x;
        zPrev = z;
        
        xPrevMod = xPrev + DELTA * zPrev;
        grad = fun.Grad(xPrevMod);
        
        x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
        z = zPrev + DELTA * (p1 * Psi(fun, xPrev) - p2 * zPrev - Hi(fun, xPrev, zPrev));
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