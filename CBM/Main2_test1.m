function Main2_test1()
    A1 = [1, 0;
          0, 1;];
    c1 = [0; 0;];
    
    A2 = [1, 0;
          0, 1;];
    c2 = [3; 4];
    
    fun1 = Func(A1, c1);
    fun2 = Func(A2, c2);
    
    spX = StartPoint(fun1, c1, c2);
    spY = StartPoint(fun2, c2, c1);
    
    CBM2(fun1, fun2, spX, spY)
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