function x = CBM_zeroing(fun, initialPoint)
    EPS = 1e-6;
    DELTA = 8;
    
    x = initialPoint;
    z = Psi(fun, x);
    
    while norm(z) > EPS
        xMod = x + DELTA * z;
        grad = fun.Grad(xMod);
        
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
        z = Psi(fun, x);
    end
end