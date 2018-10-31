function Test_2set_1()
    A_1 = [1/4, 0;
          0, 1;];
    c_1 = [3; 2;];
    
    A_2 = [1/9, 0;
          0, 1;];
    c_2 = [-4; 3];
    
    fun_1 = Func(A_1, c_1);
    fun_2 = Func(A_2, c_2);
    
    sp_X = StartPoint(fun_1, c_1, c_2);
    sp_Y = StartPoint(fun_2, c_2, c_1);
    
    tic;
    CBM_2set(fun_1, fun_2, sp_X, sp_Y);
    toc;
    
    tic;
    PM_2set(fun_1, fun_2);
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