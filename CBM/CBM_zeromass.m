function res = CBM_zeromass(fun, startingPoint)
    
end

function res = Psi(fun, x)
    gr = fun.Grad(x);
    res = (gr*(sum(x.*gr))/(norm(gr)^2) - x)/(norm(x)^3);
end

function res = Hi(fun, x, z)
    gr = fun.Grad(x);
    res = gr*(sum((fun.H()*z).*z))/(norm(gr)^2);
end
