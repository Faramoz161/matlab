function result = CBM_zeroMass_Newton(fun, initialPoint)
    % Charged balls method with zero mass
    
    EPS = 1e-6;
    DELTA = 0.5;
    
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
    gr = fun.Grad(x);
    result = p1 * fun.Psi(x) - p2 * z - gr * z.' * fun.H() * z / norm(gr)^2;
end

function result = gradF(fun, x, z)
    p2 = 1;
    gr = fun.Grad(x);
    result = - p2 * eye(length(z)) - 2 * gr * z.' * fun.H() / norm(gr)^2;
end
