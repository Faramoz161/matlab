function result = CBM_zeroMass_optimal(fun, initialPoint) % Charged balls method with zero mass and optimal delta
    EPS = 1e-6;
    
    x = initialPoint;
    z = fun.Psi(x);
    
    while norm(z) > EPS
        xPrevMod = x + OptimalDelta(fun, x, z) * z;

        grad = fun.Grad(xPrevMod);
        x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
        
        z = fun.Psi(x);
    end

    result = x;
end

function result = OptimalDelta(fun, x, z)
    delta = 50;
    lambda = 0.3;
    
    psi = norm(fun.Psi(x));
    newPsi = norm(fun.Psi(x + delta * z));
    
    while psi <= newPsi
        delta = delta * lambda;
        newPsi = norm(fun.Psi(x + delta * z));
    end
    
    result = delta;
end