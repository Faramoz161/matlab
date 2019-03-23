function result = BM_2set(fun_1, fun_2)
    % Balls method for 2 set
    
    EPS = 1e-6;
    
    A1 = fun_1.A;
    A2 = fun_2.A;
    
    normA1 = norm(A1);
    normA2 = norm(A2);
    
    c1 = fun_1.c;
    c2 = fun_2.c;
    
    while 1
        c21 = c2 - c1;
        
        ax = c21' * A1 * c21;
        ay = c21' * A2 * c21;
        
        bx = c21' * fun_1.Grad(c1);
        by = c21' * fun_2.Grad(c1);
        
        cx = fun_1.Val(c1);
        cy = fun_2.Val(c1);
        
        Dx = bx^2 - 4 * ax * cx;
        Dy = by^2 - 4 * ay * cy;
        
        t1 = (-bx + sqrt(Dx)) / (2 * ax);
        t2 = (-by - sqrt(Dy)) / (2 * ay);
        
        if t2 <= t1
            x = 0;
            y = 0;
            break;
        end
        
        x = c1 + t1 * c21;
        y = c1 + t2 * c21;
        
        temp = (x - y)' / norm(x - y);
        
        gradx = fun_1.Grad(x);
        grady = fun_2.Grad(y);
        
        tetta1 = acos(-temp * gradx / norm(gradx));
        tetta2 = acos( temp * grady / norm(grady));
        
        if (tetta1 < EPS && tetta2 < EPS)
            break;
        end
        
        c1 = x - gradx / normA1;
        c2 = y - grady / normA2;
    end
    
    result = norm(x - y);
end