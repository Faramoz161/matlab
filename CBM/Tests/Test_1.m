function Test_1()
    A = [1, 0, 0;
         0, 5, 0;
         0, 0, 5;];
    c = [11; -9; 13;];
    fun = Func(A, c);

    stz = StartPointZero(fun);
    
    tic;
    BM(fun, stz);
    toc;

    st = StartPoint(fun, c, zeros(length(c), 1));

    tic;
    CBM(fun, st);
    toc;
    
    tic;
    CBM_constStep(fun, st);
    toc;
    
    tic;
%     CBM_zeroing(fun, st);
    toc;

    tic;
    fun.Val(CBM_zeroMass(fun, st));
    toc;
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

function result = StartPointZero(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    result = t * fun.c;
end