% Charged ball method with zeroing of velocity ("optimal" step)
function result = CBM_zeroMass_optimal(fun, initialPoint) 

    EPS = 1e-6;
    
    x = initialPoint;
    z = Psi(fun, x);
    
    while norm(z) > EPS
        xMod = x + OptimalDelta(fun, x, z) * z;
        grad = fun.Grad(xMod);
        
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
        z = Psi(fun, x);
    end

    result = x;
end

function result = OptimalDelta(fun, x, z)
    delta = 50;
    lambda = 0.4;
    
    psi = norm(Psi(fun, x));
    newPsi = norm(Psi(fun, x + delta * z));
    
    while psi <= newPsi
        delta = delta * lambda;
        newPsi = norm(Psi(fun, x + delta * z));
    end
    
    result = delta;
end