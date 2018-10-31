function Graph()
    global A;
    A = [1/4 0;
           0 1/4;];
    global c;
    c = [3; 3;];
    
    initialpoint = [2; 3;];
    
    x = CBM(initialpoint);
    
    disp(x);
    disp(f(x(:,length(x))));

    filename = 'CBM.gif';
    figure(1);
    
    draw_graph();
    
    hold on;
    plot(initialpoint(1), initialpoint(2), '*');
    text(initialpoint(1) + 0.1, initialpoint(2), 'initial point');
    hold off;

    drawnow;
    
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
    
    % Add delay of 3 frames
    for i = 1 : 3
      imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
    end
    
    for i = 1 : length(x)
        draw_graph();
        
        hold on;
        plot(x(1, i), x(2, i), '*');
        
        drawnow;
        
        frame = getframe(1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im,256);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
    end
    
    text(1.6858, 1.6858, 'DONE!');
    
    drawnow;
    
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    
    % Add delay of 3 frames
    for i = 1 : 3
      imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
    end
end

function draw_graph()
    ax = [-1 4 -1 4];
    clf;
    hold on;
    
    fimplicit(@(x,y) f([x; y]), '-r', ax);
    
    fimplicit(@(x,y) x, '-k', ax);
    fimplicit(@(x,y) y, '-k', ax);

    plot(3, 3, '.m');
    plot(1.5858, 1.5858, '*c');
    
    hold off;
end

function result = f(x)
    global A;
    global c;
    result = (x - c)' * A * (x - c) - 1;
end

function result = grad(x)
    global A;
    global c;
    result = 2 * A * (x - c);
end

function result = H()
    global A;
    result = 2 * A;
end

function result = CBM(initialpoint)
    EPS = 1e-6;
    DELTA = 1;
    p1 = 10;
    p2 = 2;
    
    x = initialpoint;
    result = x;
    
    while norm(Psi(x)) > EPS
        z = p1 / p2 * Psi(x);
        xPrevMod = x + DELTA * z;
        
        gr = grad(xPrevMod);
        x = xPrevMod - gr * f(xPrevMod) / norm(gr)^2;
        
        result = [result, x];
    end
end

function result = Psi(x)
    gr = grad(x);
    result = (gr * x' * gr / norm(gr)^2 - x) / (norm(x)^3);
end