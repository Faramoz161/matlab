function res = CBM2(fun1, fun2, startX, startY)
    EPS = 1E-6;
    DELTA = 0.1;
    p1 = 10;
    p2 = 1;
    
    grad1 = fun1.Grad(startX);
    grad2 = fun2.Grad(startY);
    
    x = startX - grad1*fun1.Val(startX)/(norm(grad1)^2);
    y = startY - grad2*fun2.Val(startY)/(norm(grad2)^2);
    
    z1 = DELTA*(p1*Psi(fun1, x, y) - Hi(fun1, x, zeros(length(x), 1)));
    z2 = DELTA*(p1*Psi(fun2, y, x) - Hi(fun2, y, zeros(length(x), 1)));
    
    while norm(Psi(fun1, x, y))^2 + norm(Psi(fun2, y, x))^2 > EPS
        xPrev = x;
        yPrev = y;
        
        z1Prev = z1;
        z2Prev = z2;
        
        xPrevMod = xPrev + DELTA*z1Prev;
        yPrevMod = yPrev + DELTA*z2Prev;
        
        grad1 = fun1.Grad(xPrevMod);
        grad2 = fun2.Grad(yPrevMod);
        
        x = xPrevMod - grad1*fun1.Val(xPrevMod)/(norm(grad1)^2);
        y = yPrevMod - grad2*fun2.Val(yPrevMod)/(norm(grad2)^2);
        
        z1 = z1Prev + DELTA*(p1*Psi(fun1, xPrev, yPrev) - p2*z1Prev - Hi(fun1, xPrev, z1Prev));
        z2 = z2Prev + DELTA*(p1*Psi(fun2, yPrev, xPrev) - p2*z2Prev - Hi(fun2, yPrev, z2Prev));
    end
    
    res = norm(x - y);
end

function res = Psi(fun, x, y)
    gr = fun.Grad(x);
    temp = x - y;
    res = (gr*(sum(temp.*gr))/(norm(gr)^2) - temp)/(norm(temp)^3);
end

function res = Hi(fun, x, z)
    gr = fun.Grad(x);
    res = gr*(sum((fun.H()*z).*z))/(norm(gr)^2);
end