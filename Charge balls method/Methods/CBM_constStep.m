function result = CBM_constStep(fun, initialPoint) 
    % Charged balls method with constant step
    
    EPS = 1e-6;
    DELTA = 0.2;
    p1 = 10;
    p2 = 1;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    
    while norm(fun.Psi(x)) > EPS
        xPrev = x;
        zPrev = z;
        
        xPrevMod = xPrev + DELTA * zPrev;
        grad = fun.Grad(xPrevMod);
        
        x = xPrevMod - grad * fun.Val(xPrevMod) / norm(grad)^2;
        z = zPrev + DELTA * (p1 * fun.Psi(xPrev) - p2 * zPrev - fun.Hi(xPrev, zPrev));
    end

    result = x;
end