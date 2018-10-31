function result = CBM_zeroMass(fun, initialPoint) % Charged balls method with zero mass
    EPS = 1e-6;
    DELTA = 1;
    p1 = 10;
    p2 = 2;
    
    x = initialPoint;
    
    while norm(Psi(fun, x)) > EPS
        z = p1 / p2 * Psi(fun, x);
        xPrevMod = x + DELTA * z;

        grad = fun.Grad(xPrevMod);
        x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
    end

    result = x;
end

function result = Psi(fun, x)
    gr = fun.Grad(x);
    result = (gr * x' * gr / norm(gr)^2 - x) / norm(x)^3;
end

