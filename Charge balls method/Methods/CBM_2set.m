function result = CBM_2set(fun_1, fun_2, initial_1, initial_2) % Charged balls method for 2 set
    EPS = 1e-8;
    DELTA = 0.1;
    p1 = 10;
    p2 = 1;
    
    grad_1 = fun_1.Grad(initial_1);
    grad_2 = fun_2.Grad(initial_2);
    
    x_1 = initial_1 - grad_1 * fun_1.Val(initial_1) / norm(grad_1)^2;
    x_2 = initial_2 - grad_2 * fun_2.Val(initial_2) / norm(grad_2)^2;
    
    z_1 = DELTA * (p1 * Psi(fun_1, x_1, x_2) - Hi(fun_1, x_1, zeros(length(x_1), 1)));
    z_2 = DELTA * (p1 * Psi(fun_2, x_2, x_1) - Hi(fun_2, x_2, zeros(length(x_2), 1)));
    
    while norm(Psi(fun_1, x_1, x_2))^2 + norm(Psi(fun_2, x_2, x_1))^2 > EPS
        x_1_Prev = x_1;
        x_2_Prev = x_2;
        
        z_1_Prev = z_1;
        z_2_Prev = z_2;
        
        x_1_PrevMod = x_1_Prev + DELTA * z_1_Prev;
        x_2_PrevMod = x_2_Prev + DELTA * z_2_Prev;
        
        grad_1 = fun_1.Grad(x_1_PrevMod);
        grad_2 = fun_2.Grad(x_2_PrevMod);
        
        x_1 = x_1_PrevMod - grad_1 * fun_1.Val(x_1_PrevMod) / norm(grad_1)^2;
        x_2 = x_2_PrevMod - grad_2 * fun_2.Val(x_2_PrevMod) / norm(grad_2)^2;
        
        z_1 = z_1_Prev + DELTA * (p1 * Psi(fun_1, x_1_Prev, x_2_Prev) - p2 * z_1_Prev - Hi(fun_1, x_1_Prev, z_1_Prev));
        z_2 = z_2_Prev + DELTA * (p1 * Psi(fun_2, x_2_Prev, x_1_Prev) - p2 * z_2_Prev - Hi(fun_2, x_2_Prev, z_2_Prev));
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