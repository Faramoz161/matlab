% Polynomial approximation (an * x^n + ... + a1 * x + a0)

deg = 3;
amountPoints = [5 7 9];
f = @(x) x - sin(x) - 0.25;
a = -1;
b = 1;

for t = 1 : length(amountPoints)
    figure(t);
    hold on;
    
    n = amountPoints(t);
    x = (linspace(a, b, n));
    y = f(x);
    
    V = zeros(deg + 1, n);
    for i = 0 : deg
        V(i + 1, :) = x.^i;
    end

    coef = (V * V') \ V * y';
    
    poly = @(x) MyPolyVal(deg, coef, x);
    
    check = linspace(a, b, (b - a) * 1000);
    plot(check,    f(check), 'r-');
    plot(check, poly(check), 'k-');
    plot(    x,           y, 'b*');
    
    fmax = @(x) abs(f(x) - poly(x));
    max_dev = max(fmax(check));
    fprintf("n = %2i   max_dev = %1.6f\n", n, max_dev);
end

clear;