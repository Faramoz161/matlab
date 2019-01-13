function result = CBM_zeroing(fun, initialPoint)
    % Charged balls method with zeroing speed
    
    EPS = 1e-6;
    DELTA = 0.25;
    p1 = 100;
    
    x = initialPoint;
    z = DELTA * (p1 * fun.Psi(x) - fun.Hi(x, zeros(length(x), 1)));
    
    while norm(fun.Psi(x)) > EPS
        xPrev = x;
        xPrevMod = xPrev + DELTA * z;
        grad = fun.Grad(xPrevMod);
        x = xPrevMod - grad * fun.Val(xPrevMod) / (norm(grad)^2);
        z = DELTA * (p1 * fun.Psi(x) - fun.Hi(xPrev, zeros(length(x), 1)));
    end
    
    result = x;
end