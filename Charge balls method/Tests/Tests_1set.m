function Tests_1set()
    time_BM = 0;
    time_CBM_constStep = 0;
    time_CBM_zeroMass = 0;
    time_CBM_zeroMass_Newton = 0;
    time_CBM_zeroMass_optimal = 0;
    time_MyMethod = 0;
    
    dim = 1000;
    A = diag(1 + 4 * rand(dim, 1));
    
    for amount = 1 : 1
        c = 2 * rand(dim, 1) - 1;
        c = c / norm(c) * 6;
        
        fun = Func(A, c);
        st = StartPointZero(fun);
        
%         tic;
%         BM(fun, st);
%         time_BM = time_BM + toc;
        
%         tic;
%         CBM_constStep(fun, st);
%         time_CBM_constStep = time_CBM_constStep + toc;
        
        tic;
        CBM_zeroMass(fun, st);
        time_CBM_zeroMass = time_CBM_zeroMass + toc;
        
        tic;
        CBM_zeroMass_Newton(fun, st);
        time_CBM_zeroMass_Newton = time_CBM_zeroMass_Newton + toc;
        
        tic;
        CBM_zeroMass_optimal(fun, st);
        time_CBM_zeroMass_optimal = time_CBM_zeroMass_optimal + toc;
        
        tic;
        MyMethod(fun, st);
        time_MyMethod = time_MyMethod + toc;
    end
    
%     fprintf("%g\n", time_BM);
%     fprintf("%g\n", time_CBM_constStep);
    fprintf("%g\n", time_CBM_zeroMass);
    fprintf("%g\n", time_CBM_zeroMass_Newton);
    fprintf("%g\n", time_CBM_zeroMass_optimal);
    fprintf("%g\n", time_MyMethod);
end

function result = StartPointZero(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    result = t * fun.c;
end