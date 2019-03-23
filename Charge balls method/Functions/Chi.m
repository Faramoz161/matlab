% Calculation chi(x, z)
function result = Chi(fun, x, z)
    grad = fun.Grad(x);
    result = grad * z' * fun.H() * z / norm(grad)^2;
end