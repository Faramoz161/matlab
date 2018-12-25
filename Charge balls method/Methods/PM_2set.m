function result = PM_2set(fun_1, fun_2)
    % Penalty method for 2 set
    
    EPS = 1e-6;
    
    n = length(fun_1.c);
    z = zeros(2 * n, 1);
    r = 1;
    
    while H(fun_1, fun_2, z) > EPS
        
        while norm(F_grad(fun_1, fun_2, z, r)) > EPS / 10
            z = z - F_H(fun_1, fun_2, z, r) \ F_grad(fun_1, fun_2, z, r);
        end

        r = r * 10;
    end
    
    x = z(1 : n);
    y = z(n + 1 : 2 * n);
    
    result = norm(x - y);
end

function result = F_grad(fun_1, fun_2, z, r)
    n = length(z) / 2;
    x = z(1 : n);
    y = z(n + 1 : 2 * n);
    
    term_1 = 2 * [(x - y);
                  (y - x);];
    term_2 = 2 * r * [fun_1.Val(x) * fun_1.Grad(x); 
                      fun_2.Val(y) * fun_2.Grad(y);];
    
    result = term_1 + term_2;
end

function result = F_H(fun_1, fun_2, z, r)
    n = length(z) / 2;
    x = z(1 : n);
    y = z(n + 1 : 2 * n);
    
    grad_1 = fun_1.Grad(x);
    grad_2 = fun_2.Grad(y);
    
    term_1 = 2 * [ eye(n), -eye(n); 
                  -eye(n),  eye(n);];
    term_2 = 2 * r * [grad_1 * grad_1' + fun_1.Val(x) * fun_1.H(), zeros(n);
                      zeros(n), grad_2 * grad_2' + fun_2.Val(y) * fun_2.H();];
    
    result = term_1 + term_2;
end

function result = H(fun_1, fun_2, z)
    n = length(z) / 2;
    x = z(1 : n);
    y = z(n + 1 : 2 * n);
    
    result = fun_1.Val(x)^2 + fun_2.Val(y)^2;
end