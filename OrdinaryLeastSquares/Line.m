% Line approximation (a1 * x + a0)

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
    
    sumX = sum(x);
    sumX2 = sum(x * x');
    sumY = sum(y);
    sumYX = sum(y * x');
    
    a0 = (sumX2 * sumY - sumYX * sumX) / (n * sumX2 - sumX^2);
    a1 = (n * sumYX - sumY * sumX) / (n * sumX2 - sumX^2);

    line = @(x) a1 * x + a0;
    
    m = linspace(a, b, (b - a) * 1000);
    plot(m, line(m), 'k-');
    plot(m,    f(m), 'r-');
    plot(x,       y, 'b*');
    
    fmax = @(x) abs(f(x) - line(x));
    max_dev = max(fmax(m));
    fprintf("n = %2i   max_dev = %1.6f\n", n, max_dev);
end

clear;