function result = CBM_2set_zm(fun_1, fun_2, initial_1, initial_2) % Charged balls method for 2 set
    EPS = 1e-6;
    DELTA = 10;
    
    x_1 = initial_1;
    x_2 = initial_2;
    
    while norm(Psi(fun_1, x_1, x_2))^2 + norm(Psi(fun_2, x_2, x_1))^2 > EPS
        z_1 = Psi(fun_1, x_1, x_2);
        z_2 = Psi(fun_2, x_2, x_1);
        
        x_1_PrevMod = x_1 + DELTA * z_1;
        x_2_PrevMod = x_2 + DELTA * z_2;
        
        grad_1 = fun_1.Grad(x_1_PrevMod);
        grad_2 = fun_2.Grad(x_2_PrevMod);
        
        x_1 = x_1_PrevMod - grad_1 * fun_1.Val(x_1_PrevMod) / norm(grad_1)^2;
        x_2 = x_2_PrevMod - grad_2 * fun_2.Val(x_2_PrevMod) / norm(grad_2)^2;
    end
    
    result = norm(x_1 - x_2);
end

function result = Psi(fun, x_1, x_2)
    gr = fun.Grad(x_1);
    temp = x_1 - x_2;
    result = (gr * temp' * gr / norm(gr)^2 - temp) / norm(temp)^3;
end

function result = Hi(fun, x, z)
    gr = fun.Grad(x);
    result = gr * z' * fun.H() * z / norm(gr)^2;
end