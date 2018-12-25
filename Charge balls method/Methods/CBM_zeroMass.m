function result = CBM_zeroMass(fun, initialPoint)
    % Charged balls method with zero mass
    
    EPS = 1e-6;
    DELTA = 10;
    
    x = initialPoint;
    z = fun.Psi(x);
    
    while norm(z) > EPS
        xMod = x + DELTA * z;
        grad = fun.Grad(xMod);
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
        z = fun.Psi(x);
        
        DELTA = DELTA * 0.99;
    end

    result = x;
end