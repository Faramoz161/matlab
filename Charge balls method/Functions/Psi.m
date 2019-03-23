% Calculation psi(x)
function result = Psi(fun, x)
    grad = fun.Grad(x);
    result = (grad * x' * grad / norm(grad)^2 - x) / norm(x)^3;
end