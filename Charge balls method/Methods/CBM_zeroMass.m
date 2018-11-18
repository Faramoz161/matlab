function result = CBM_zeroMass(fun, initialPoint) % Charged balls method with zero mass
    EPS = 1e-6;
    DELTA = 10;
    
    x = initialPoint;
    
    while norm(fun.Psi(x)) > EPS
        z = fun.Psi(x);
        xPrevMod = x + DELTA * z;

        grad = fun.Grad(xPrevMod);
        x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
    end

    result = x;
end