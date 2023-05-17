function P = fullify(o, roundparam)
    % Ax=b, lb <= x <= ub
    A = o.A;
    b = o.b;
    lb = o.barrier.lb;
    ub = o.barrier.ub;
    x = o.center;
    
    % shift the center to 0: the problem becomes Ax = 0, lb <= x <= ub
    lb = lb - x;
    ub = ub - x;
    
    % rescale the polytope
    s = min(-lb, ub);
    A = A .* s';
    lb = lb ./ s;
    ub = ub ./ s;
    
    P.s = s;
    P.A_eq = A; P.b_eq = zeros(size(A, 1), 1);
    P.center = x;
    P.b = [ub; -lb];
    
    tic;
    P = newPreprocess(P, roundparam);
    P.T = o.T; P.y = o.y;
    roundTime = toc;
    P.roundTime = roundTime;
end