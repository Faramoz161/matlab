% Spline 3,2

amountPoints = [3 5 8 10];
f = @(x) x - sin(x) - 0.25;
a = -1;
b = 1;

for t = 1 : length(amountPoints)
    figure(t);
    hold on;
    
    n = amountPoints(t);
    x = linspace(a, b, n);
    y = f(x);
    
    h = x(2) - x(1);
    
    H = 4 * h * eye(n - 2);
    for i = 1 : n - 3
        H(i, i + 1) = h;
        H(i + 1, i) = h;
    end
    
    gamma = zeros(n - 2, 1);
    for i = 2 : n - 1
        gamma(i - 1) = 6 / h * (y(i + 1) - 2 * y(i) + y(i - 1));
    end
    
    y_2 = [sin(a); H \ gamma; sin(b);];
    
    y_1 = zeros(n - 1, 1);
    for i = 1 : n - 1
        y_1(i) = (y(i + 1) - y(i)) / h - y_2(i + 1) * h / 6 - y_2(i) * h / 3;
    end
    
    max_dev = 0;
    for i = 1 : n - 1
        spline = @(z) y(i) + y_1(i) * (z - x(i)) + y_2(i) * (z - x(i)).^2 / 2 + (y_2(i + 1) - y_2(i)) * (z - x(i)).^3 / 6 / h;
        check = linspace(x(i), x(i + 1));
        plot(check,      f(check), 'r-');
        plot(check, spline(check), 'k-');

        fmax = @(x) abs(f(x) - spline(x));
        max_dev = max(max(fmax(check)), max_dev);
    end
    
    fprintf("n = %2i   max_dev = %1.6f\n", n, max_dev);
    plot(x, y, 'b*');
end

clear;