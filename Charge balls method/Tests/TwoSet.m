function TwoSet()
    time_CBM = 0;
    time_PM = 0;
    
    dim = 3;
    A_1 = diag(1 + 4 * rand(dim, 1));
    A_2 = diag(1 + 4 * rand(dim, 1));
    
    for amount = 1 : 1000
        c = 2 * (rand(dim, 1) - 0.5);
        distance = 5.1;
        c = c / norm(c) * distance;
        
        fun_1 = Func(A_1, c);
        fun_2 = Func(A_2, -c);

        sp_X = StartPointZero(fun_1);
        sp_Y = StartPointZero(fun_2);

        tic;
        CBM_2set_zm(fun_1, fun_2, sp_X, sp_Y);
        time_CBM = time_CBM + toc;

        tic;
        PM_2set(fun_1, fun_2);
        time_PM = time_PM + toc;
    end
    
    fprintf("CBM_2set_zm: %g\n", time_CBM / amount);
    fprintf("PM_2set:     %g\n", time_PM / amount);
end

function result = StartPointZero(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    result = t * fun.c;
end