function Main()
    A = [1, 0, 0;
         0, 5, 0;
         0, 0, 5;];
%     c = [11; -9; 13;];
    
    s = 0;
    e = pi/2;
    step = pi/4;
    distance = 3;
    
    time_BM = 0;
    time_CBM = 0;
    
    for psi = s : step : e
        for fi = s : step : e
            c = distance*[cos(fi)*cos(psi); sin(fi)*cos(psi); sin(psi);];
%             c = [0; 0; 3;];
            fun = Func(A, c);
            
            stz = StartPointZero(fun);
            tic;
            BM(fun, stz);
            time_BM = time_BM + toc;

            st = StartPoint(fun, c, zeros(length(c), 1));
            tic;
            CBM(fun, st);
            time_CBM = time_CBM + toc;
            
            if psi == e
                break;
            end
        end
    end
    
    disp(time_BM);
    disp(time_CBM);
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