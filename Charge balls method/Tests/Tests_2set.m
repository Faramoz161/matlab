function Tests_2set()
    dimension = 3;
    time_CBM = 0;
    time_PM = 0;
    
    A_1 = eye(dimension);
    for i = 1 : dimension
        A_1(i,i) = 1 + 4 * rand();
    end
    
    A_2 = eye(dimension);
    for i = 1 : dimension
        A_2(i,i) = 1 + 4 * rand();
    end
    
    for amount = 1 : 1
        c = 2 * (rand(dimension, 1) - 0.5);
        distance = 5.1 + 4 * rand();
        c = c / norm(c) * distance;
        
        fun_1 = Func(A_1, c);
        fun_2 = Func(A_2, -c);

        sp_X = StartPoint(fun_1, c, zeros(dimension, 1));
        sp_Y = StartPoint(fun_2, -c, zeros(dimension, 1));

        tic;
        CBM_2set_zm(fun_1, fun_2, sp_X, sp_Y)
        time_CBM = time_CBM + toc;

        tic;
        PM_2set(fun_1, fun_2)
        time_PM = time_PM + toc;
    end
    
    fprintf('%g\n', time_CBM / amount);
    fprintf('%g\n', time_PM / amount);
end

function result = StartPoint(fun, left, right)
    while norm(left - right) > 1e-8
        mid = (left + right) / 2;
        if (fun.Val(mid) < 0)
            left = mid;
        else
            right = mid;
        end
    end
    
    result = left;
end