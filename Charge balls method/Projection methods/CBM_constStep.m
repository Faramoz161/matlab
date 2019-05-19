function x = CBM_constStep(fun, initialPoint)
    EPS = 1e-6;
    DELTA = 0.7;
    p1 = 30;
    p2 = 1;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    psi = Psi(fun, x);
    
    while norm(psi) > EPS
        z = z + DELTA * (p1 * psi - p2 * z - Chi(fun, x, z));
        
        xMod = x + DELTA * z;
        grad = fun.Grad(xMod);
        
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
        psi = Psi(fun, x);
    end
end