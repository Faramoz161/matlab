function res = CBM(fun, startingPoint)
    EPS = 1E-6;
    EPS1 = 1E-8;
    EPS2 = 1E-2;
    DELTA = 0.1;
    p1 = 100;
    p2 = 2;
    
    grad = fun.Grad(startingPoint);
    x = startingPoint - grad*fun.Val(startingPoint)/(norm(grad)^2);
    z = DELTA*(p1*Psi(fun, x) - Hi(fun, x, zeros(length(x), 1)));
    
    while norm(Psi(fun, x)) > EPS
        xPrev = x;
        zPrev = z;

        %while 1
            xPrevMod = xPrev + DELTA*zPrev;
            grad = fun.Grad(xPrevMod);
            x = xPrevMod - grad*fun.Val(xPrevMod)/(norm(grad)^2);
            
        %    condition = norm(x - xPrevMod)/norm(xPrevMod);
        %    if condition < EPS1
        %        DELTA = DELTA*2;
        %    elseif condition > EPS2
        %        DELTA = DELTA/2;
        %    else
        %        break;
        %    end
        %end
        
        z = zPrev + DELTA*(p1*Psi(fun, xPrev) - p2*zPrev - Hi(fun, xPrev, zPrev));
    end
    
    res = norm(x);
end

function res = Psi(fun, x)
    gr = fun.Grad(x);
    res = (gr*(sum(x.*gr))/(norm(gr)^2) - x)/(norm(x)^3);
end

function res = Hi(fun, x, z)
    gr = fun.Grad(x);
    res = gr*(sum((fun.H()*z).*z))/(norm(gr)^2);
end