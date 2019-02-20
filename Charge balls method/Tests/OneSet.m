function OneSet()
    time_BM                   = 0;
    time_CBM_constStep        = 0;
    time_CBM_zeroing          = 0;
    time_CBM_zeroMass         = 0;
    time_CBM_zeroMass_Newton  = 0;
    time_CBM_zeroMass_optimal = 0;
    time_PM                   = 0;
    
    dim = 3;
    A = diag(1 + 4 * rand(dim, 1));
    
    for amount = 1 : 1000
        c = 2 * rand(dim, 1) - 1;
        distance = 5.5;
        c = c / norm(c) * distance;

        fun = Func(A, c);
        startPoint = StartPointZero(fun);

        tic;
        BM(fun, startPoint);
        time_BM = time_BM + toc;
        
        tic;
        CBM_constStep(fun, startPoint);
        time_CBM_constStep = time_CBM_constStep + toc;

        tic;
        CBM_zeroing(fun, startPoint);
        time_CBM_zeroing = time_CBM_zeroing + toc;
        
        tic;
        CBM_zeroMass(fun, startPoint);
        time_CBM_zeroMass = time_CBM_zeroMass + toc;

        tic;
        CBM_zeroMass_Newton(fun, startPoint);
        time_CBM_zeroMass_Newton = time_CBM_zeroMass_Newton + toc;
        
        tic;
        CBM_zeroMass_optimal(fun, startPoint);
        time_CBM_zeroMass_optimal = time_CBM_zeroMass_optimal + toc;
        
        tic;
        PM_1set(fun, startPoint);
        time_PM = time_PM + toc;
    end
    
    fprintf("BM:                   %g\n", time_BM / amount);
    fprintf("CBM_constStep:        %g\n", time_CBM_constStep / amount);
    fprintf("CBM_zeroing:          %g\n", time_CBM_zeroing / amount);
    fprintf("CBM_zeroMass:         %g\n", time_CBM_zeroMass / amount);
    fprintf("CBM_zeroMass_Newton:  %g\n", time_CBM_zeroMass_Newton / amount);
    fprintf("CBM_zeroMass_optimal: %g\n", time_CBM_zeroMass_optimal / amount);
    fprintf("PM:                   %g\n", time_PM / amount);
end

function result = StartPointZero(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    result = t * fun.c;
end