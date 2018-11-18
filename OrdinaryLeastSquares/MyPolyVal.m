function result = MyPolyVal(deg, coef, x)
    result = 0;
    for i = 0 : deg
        result = result + coef(i + 1) * x.^i;
    end
end

