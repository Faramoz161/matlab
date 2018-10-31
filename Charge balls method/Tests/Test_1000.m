function Test_1000()
    dimension = 1000;
    distance = 5.5;
    time_BM = 0;
    time_CBM = 0;
    time_CBM_constStep = 0;
%     time_CBM_zeroing = 0;
    time_CBM_zeroMass = 0;
    
    A = eye(dimension);
    for i = 1 : dimension
        A(i,i) = 1 + 4 * rand();
    end
    
    for amount = 1 : 1
        c = 2 * (rand(dimension, 1) - 0.5);
        c = c / norm(c) * distance;
        
        fun = Func(A, c);
        
        st = StartPointZero(fun);
        
        tic;
        BM(fun, st);
        time_BM = time_BM + toc;
        
%         tic;
%         CBM(fun, st);
%         time_CBM = time_CBM + toc;
        
%         tic;
%         CBM_constStep(fun, st);
%         time_CBM_constStep = time_CBM_constStep + toc;
        
%         tic;
%         CBM_zeroing(fun, st);
%         time_CBM_zeroing = time_CBM_zeroing + toc;
        
        tic;
        CBM_zeroMass(fun, st);
        time_CBM_zeroMass = time_CBM_zeroMass + toc;
    end
    
    fprintf('\n%d\n', time_BM);
%     fprintf('%d\n', time_CBM);
%     fprintf('%d\n', time_CBM_constStep);
%     fprintf('%d\n', time_CBM_zeroing);
    fprintf('%d\n', time_CBM_zeroMass);
end

function result = StartPointZero(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    result = t * fun.c;
end