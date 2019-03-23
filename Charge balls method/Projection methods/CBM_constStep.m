% Charged balls method with constant step
function result = CBM_constStep(fun, initialPoint)

    EPS = 1e-6;
    DELTA = 0.3;
    p1 = 30;
    p2 = 1;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    
    while norm(Psi(fun, x)) > EPS
        xMod = x + DELTA * z;
        grad = fun.Grad(xMod);
        
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
        z = z + DELTA * (p1 * Psi(fun, x) - p2 * z - Chi(fun, x, z));
    end

    result = x;
end