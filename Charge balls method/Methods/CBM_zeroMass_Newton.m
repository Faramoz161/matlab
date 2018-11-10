function result = CBM_zeroMass_Newton(fun, initialPoint)
    EPS = 1e-6;
    DELTA = 1;
    
    x = initialPoint;
    z = zeros(length(x), 1);
    
    while norm(Psi(fun, x)) > EPS
        z = z - gradF(fun, x, z) \ F(fun, x ,z); %NewtonStep
        xPrev = x;
        
        xPrevMod = xPrev + DELTA*z;
        grad = fun.Grad(xPrevMod);
        
        x = xPrevMod - grad*fun.Val(xPrevMod)/(norm(grad)^2);
    end
    
    result = x;
end

function res = Psi(fun, x)
    gr = fun.Grad(x);
    res = (gr*(sum(x.*gr))/(norm(gr)^2) - x)/(norm(x)^3);
end

function res = F(fun, x, z)
    p1 = 10;
    p2 = 1;
    gr = fun.Grad(x);
    res = p1 * Psi(fun, x) - p2 * z - gr * z.' * fun.H() * z / (norm(gr)^2);
end

function res = gradF(fun, x, z)
    p2 = 1;
    gr = fun.Grad(x);
    res = - p2 * eye(length(z)) - 2 * gr * z.' * fun.H() / (norm(gr)^2);
end
