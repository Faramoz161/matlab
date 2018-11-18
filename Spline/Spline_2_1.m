% Spline 2,1

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
    
    max_dev = 0;
    a_prev_2 = 0;
    a_prev_1 = 1 - cos(a);
    
    for i = 1 : n - 1
        X = [    x(i)^2      x(i)  1; 
             x(i + 1)^2  x(i + 1)  1;
               2 * x(i)         1  0;];
        Y = [    y(i);
             y(i + 1);
             2 * x(i) * a_prev_2 + a_prev_1;];
        A = X \ Y;
        
        a_prev_2 = A(1);
        a_prev_1 = A(2);
        
        spline = @(x) A(1) * x.^2 + A(2) * x + A(3);
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