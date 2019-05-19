function x = CBM_zeroing_frag(fun, initialPoint) 
    EPS = 1e-6;
    
    x = initialPoint;
    z = Psi(fun, x);
    
    while norm(z) > EPS
        xMod = x + OptimalDelta(fun, x, z) * z;
        grad = fun.Grad(xMod);
        
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
        z = Psi(fun, x);
    end
end

function delta = OptimalDelta(fun, x, z)
    delta = 50;
    lambda = 0.4;
    
    psi = norm(Psi(fun, x));
    newPsi = norm(Psi(fun, x + delta * z));
    
    while psi <= newPsi
        delta = delta * lambda;
        newPsi = norm(Psi(fun, x + delta * z));
    end
end