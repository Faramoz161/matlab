function PointProjection()
    % ------------------------------------
    % Testing methods for point projection
    % ------------------------------------
    
    % Set dimensions and number of tests
    dimensions = [2; 3; 10; 100; 500; 1000];
    testsCount = 100;
    
    % Distance to the center of the ellipse
    distance = 5;
    
    % Time setting
    zero = zeros(length(dimensions), 1);
    time_BM               = zero;
    time_CBM              = zero;
    time_CBM_zeroing      = zero;
    time_CBM_zeroMass     = zero;
    time_CBM_zeroing_frag = zero;
    time_PM               = zero;
    
    for i = 1 : length(dimensions)
        dim = dimensions(i);
        
        for j = 1 : testsCount
            % Matrix A of the ellipse
            A = diag(2 + 3 * rand(dim, 1));
            
            % Center of thera ellipse
            c = normrnd(0, 1, dim, 1);
            c = c / norm(c) * distance;
            
            % Set the ellipse
            fun = Func(A, c);
            
            % Calculate the starting point
            startpoint = StartingPoint(fun);
    
            % Charged ball method with constant step
            tic;
            CBM_constStep(fun, startpoint);
            time_CBM(i) = time_CBM(i) + toc;
            
            % Charged ball method with zero mass (using Newton's method)
            tic;
            CBM_zeroMass(fun, startpoint);
            time_CBM_zeroMass(i) = time_CBM_zeroMass(i) + toc;
            
            % Charged ball method with zeroing of velocity
            tic;
            CBM_zeroing(fun, startpoint);
            time_CBM_zeroing(i) = time_CBM_zeroing(i) + toc;

            % Charged ball method with zeroing of velocity (fragmentation step)
            tic;
            CBM_zeroing_frag(fun, startpoint);
            time_CBM_zeroing_frag(i) = time_CBM_zeroing_frag(i) + toc;
            
            % Penalty method
            tic;
            PM_projection(fun, startpoint);
            time_PM(i) = time_PM(i) + toc;
            
            % Ball method
            tic;
            BM(fun, startpoint);
            time_BM(i) = time_BM(i) + toc;
        end
    end
    
    % Average time
    time_CBM = time_CBM / testsCount;
    time_CBM_zeroMass = time_CBM_zeroMass / testsCount;
    time_CBM_zeroing = time_CBM_zeroing / testsCount;
    time_CBM_zeroing_frag = time_CBM_zeroing_frag / testsCount;
    time_PM = time_PM / testsCount;
    time_BM = time_BM / testsCount;
        
    % Create table and display
    output = table(dimensions, time_CBM, time_CBM_zeroMass, ...
                   time_CBM_zeroing, time_CBM_zeroing_frag, ...
                   time_PM, time_BM);
    disp(output);

    f = figure;

    hold on;
    grid on;

    loglog(dimensions, time_CBM, '-s');
    loglog(dimensions, time_CBM_zeroMass, '-s');
    loglog(dimensions, time_CBM_zeroing, '-s');
    loglog(dimensions, time_CBM_zeroing_frag, '-s');
    loglog(dimensions, time_PM, '-s');
    loglog(dimensions, time_BM, '-s');
    
    legend('CBM', 'CBM_m', 'CBM_v', 'CBM_{vf}', 'PM', 'BM');
    
    f.CurrentAxes.XScale = 'log';
    f.CurrentAxes.YScale = 'log';
end

function result = StartingPoint(fun)
    t = 1 - 1 / sqrt(fun.c' * fun.A * fun.c);
    result = t * fun.c;
end