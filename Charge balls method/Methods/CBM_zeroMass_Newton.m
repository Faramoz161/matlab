function result = CBM_zeroMass_Newton(fun, initialPoint)
    % Charged balls method with zero mass
    
    EPS = 1e-6;
    DELTA = 1;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    
    while norm(fun.Psi(x)) > EPS
        z = z - gradF(fun, x, z) \ F(fun, x ,z); % Newton step
        
        xMod = x + DELTA * z;
        grad = fun.Grad(xMod);
        x = xMod - grad * fun.Val(xMod) / norm(grad)^2;
    end
    
    result = x;
end

function result = F(fun, x, z)
    p1 = 10;
    p2 = 1;
    result = p1 * fun.Psi(x) - p2 * z - fun.Hi(x, z);
end

function result = gradF(fun, x, z)
    p2 = 1;
    gr = fun.Grad(x);
    result = - p2 * eye(length(z)) - 2 * gr * z.' * fun.H() / norm(gr)^2;
end
