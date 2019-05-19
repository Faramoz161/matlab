function x = CBM_zeroMass(fun, initialPoint)
    EPS = 1e-6;
    DELTA = 0.8;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    
    while norm(Psi(fun, x)) > EPS
        z = z - gradF(fun, x, z) \ F(fun, x ,z);
        
        xMod = x + DELTA * z;
        grad = fun.Grad(xMod);
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
    end
end

function result = F(fun, x, z)
    p1 = 10;
    p2 = 1;
    result = p1 * Psi(fun, x) - p2 * z - Chi(fun, x, z);
end

function result = gradF(fun, x, z)
    p2 = 1;
    gr = fun.Grad(x);
    result = - p2 * eye(length(z)) - 2 * gr * z.' * fun.H() / norm(gr)^2;
end
