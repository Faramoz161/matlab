function PointProjection()
    % ------------------------------------
    % Testing methods for point projection
    % ------------------------------------
    
    % Set dimensions and number of tests
    dimensions = [2, 3, 5, 7, 10]';
    testsCount = 1000;
    
    % Distance to the center of the ellipse
    distance = 5;
    
    % Time setting
    zero = zeros(length(dimensions), 1);
    time_BM                      = zero;
    time_CBM_constStep           = zero;
    time_CBM_zeroing             = zero;
    time_CBM_zeroMass_Newton     = zero;
    time_CBM_zeroing_optimalStep = zero;
    time_PM                      = zero;
    
    for i = 1 : length(dimensions)
        dim = dimensions(i);
        
        for j = 1 : testsCount
            % Matrix A of the ellipse
            A = diag(2 + 3 * rand(dim, 1));
            
            % Center of the ellipse
            c = rand(dim, 1) - 0.5;
            c = c / norm(c) * distance;
            
            % Set the ellipse
            fun = Func(A, c);
            
            % Calculate the starting point
            startpnt = StartingPoint(fun);
            
            % Ball method
            tic;
            BM(fun, startpnt);
            time_BM(i) = time_BM(i) + toc;
    
            % Charged ball method with constant step
            tic;
            CBM_constStep(fun, startpnt);
            time_CBM_constStep(i) = time_CBM_constStep(i) + toc;
            
            % Charged ball method with zeroing of velocity
            tic;
            CBM_zeroing(fun, startpnt);
            time_CBM_zeroing(i) = time_CBM_zeroing(i) + toc;

            % Charged ball method with zero mass (using Newton's method)
            tic;
            CBM_zeroMass_Newton(fun, startpnt);
            time_CBM_zeroMass_Newton(i) = time_CBM_zeroMass_Newton(i) + toc;
            
            % Charged ball method with zeroing of velocity ("optimal" step)
            tic;
            CBM_zeroMass_optimal(fun, startpnt);
            time_CBM_zeroing_optimalStep(i) = time_CBM_zeroing_optimalStep(i) + toc;
            
            % Penalty method
            tic;
            PM_projection(fun, startpnt);
            time_PM(i) = time_PM(i) + toc;
        end
    end
    
    % Average time
    time_CBM_constStep = time_CBM_constStep / testsCount;
    time_CBM_zeroing = time_CBM_zeroing / testsCount;
    time_CBM_zeroing_optimalStep = time_CBM_zeroing_optimalStep / testsCount;
    time_CBM_zeroMass_Newton = time_CBM_zeroMass_Newton / testsCount;
    time_BM = time_BM / testsCount;
    time_PM = time_PM / testsCount;
    
    % Create table and display
    output = table(dimensions, time_CBM_constStep, time_CBM_zeroing, ...
                   time_CBM_zeroing_optimalStep, time_CBM_zeroMass_Newton, ...
                   time_BM, time_PM);
    disp(output);    
end

function sp = StartingPoint(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    sp = t * fun.c;
end