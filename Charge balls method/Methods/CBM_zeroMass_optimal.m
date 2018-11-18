function result = CBM_zeroMass_optimal(fun, initialPoint) % Charged balls method with zero mass
    EPS = 1e-6;
    
    x = initialPoint;
    z = Psi(fun, x);
    
    while norm(z) > EPS
        xPrevMod = x + OptimalDelta(fun, x, z) * z;

        grad = fun.Grad(xPrevMod);
        x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
        
        z = Psi(fun, x);
    end

    result = x;
end

function result = Psi(fun, x)
    gr = fun.Grad(x);
    result = (gr * x' * gr / norm(gr)^2 - x) / norm(x)^3;
end

function result = OptimalDelta(fun, x, z) %Golden ratio
    LEFT_CONST = (3 - sqrt(5)) / 2;
    RIGHT_CONST = (sqrt(5) - 1) / 2;
    
    b = 65;
    c = LEFT_CONST * b;
    d = RIGHT_CONST * b;
    
    psiC = Psi(fun, x + c * z);
    psiD = Psi(fun, x + d * z);
    
    while b > 10
        if norm(psiC) <= norm(psiD)
            b = d;            
            d = c;
            c = LEFT_CONST * b;
            
            psiD = psiC;
            psiC = Psi(fun, x + c * z);
        else
            result = d;
            return;
        end
    end
    result = b;
end