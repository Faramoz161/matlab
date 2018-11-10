function result = MyMethod(fun, initialPoint)
    EPS = 1e-6;
    
    x = initialPoint;
    z = Psi(fun, x);
    
    while norm(z) > EPS
        xZero = x + OptimalDelta(fun, x, z) * z;
    
        xZero = xZero / norm(xZero) * 6;
        x = PointZero(fun, xZero, zeros(length(xZero), 1));
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
    
    a = 0;
    b = 100;
    c = LEFT_CONST * (b - a) + a;
    d = RIGHT_CONST * (b - a) + a;
    
    psiC = Psi(fun, x + c * z);
    psiD = Psi(fun, x + d * z);
    
    while b - a > 2e-3
        if norm(psiC) <= norm(psiD)
            b = d;            
            d = c;
            c = LEFT_CONST * (b - a) + a;
            
            psiD = psiC;
            psiC = Psi(fun, x + c * z);
        else
            a = c;
            c = d;
            d = RIGHT_CONST * (b - a) + a;
            
            psiC = psiD;
            psiD = Psi(fun, x + d * z);
        end
    end
    
    result = (b + a) / 2;
end

function result = PointZero(fun, left, right)
    while norm(left - right) > 1e-3
        mid = (left + right) / 2;
        if (fun.Val(mid) < 0)
            left = mid;
        else
            right = mid;
        end
    end
    
    result = left;
end