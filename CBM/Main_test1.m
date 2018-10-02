function Main_test1()
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
end

function res = StartPoint(fun, left, right)
    while (norm(left - right) > 1E-8)
        mid = (left + right)/2;
        if (fun.Val(mid) < 0)
            left = mid;
        else
            right = mid;
        end
    end
    
    res = left;
end

function res = StartPointZero(fun)
    t = 1 - 1/(sqrt(sum(fun.c .* (fun.A*fun.c))));
    res = t*fun.c;
end