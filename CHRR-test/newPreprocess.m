function [P] = newPreprocess(P, roundparam)
%PREPROCESS
% INPUTS:
% P ... polytope object
% .A ...
% .b ... specifying Ax<=b
% .A_eq ...
% .b_eq ... specifying A_eq x = b_eq
% roundparam    0: no rounding
%               1 (default): round using max volume ellipsoid

%restrict to the degenerate subspace, if it exists
N = null(full(P.A_eq));

%find a point in the null space to define the shift
p_shift = P.center;
P.A = [N; -N];
N_total = P.s .* N;

%remove zero rows
row_norms = sqrt(sum(P.A.^2,2));
P.A = P.A(row_norms>1e-6,:);
P.b = P.b(row_norms>1e-6,:);
dim = size(P.A,2);
fprintf('Now in %d dimensions after restricting\n', dim);

%scale the matrix to help with numerical precision
if exist('gmscale')==2
    fprintf('Preconditioning A with gmscale\n');
    [cs,rs] = gmscale(P.A,0,0.99);
    P.A = diag(1./rs)*P.A*diag(1./cs);
    P.b = diag(1./rs)*P.b;
    N_total = N_total*diag(1./cs);
end

row_norms = sqrt(sum(P.A.^2,2));
P.A = diag(1./row_norms)*P.A;
P.b = diag(1./row_norms)*P.b;

if roundparam==1
    fprintf('Rounding...\n');
    T = eye(dim);

    if dim<200
        x0 = linprog(zeros(size(P.A, 2), 1), P.A, P.b);
        
        [T_shift, Tmve, converged] = mve_run_cobra(P.A, P.b, x0, 1e-8);
        [P, N_total, p_shift, T] = shiftPolytope(P, N_total, p_shift, T, Tmve, T_shift);
        if converged~=1
            fprintf('There was a problem with finding the maximum volume ellipsoid.\n');
        end
    else
        %if dimension is large, be a little more careful
        %the below loop looks silly for an interior point method, but is actually
        %quite important for numerical stability. while normally you'd only call an optimization routine
        %once, we call it iteratively--where in each iteration, we map a large
        %ellipsoid to the unit ball. the idea is that the "iterative roundings"
        %will make each subsequent problem easier to solve numerically.
        max_its = 20;
        its = 0;
        reg = 1e-3;
        Tmve = eye(dim);
        converged = 0;
        while (max(eig(Tmve))>6*min(eig(Tmve)) && converged~=1) || reg>1e-6 || converged==2
            its = its+1;
            x0 = linprog(zeros(size(P.A, 2), 1), P.A, P.b);

            reg = max(reg/10,1e-10);
            [T_shift, Tmve, converged] = mve_run_cobra(P.A, P.b, x0, reg);
            [P, N_total, p_shift, T] = shiftPolytope(P, N_total, p_shift, T, Tmve, T_shift);
            row_norms = sqrt(sum(P.A.^2,2));
            P.A = diag(1./row_norms)*P.A;
            P.b = diag(1./row_norms)*P.b;
            if its==max_its
                break;
            end

            fprintf('Iteration %d: reg=%.1e, ellipsoid vol=%.1e, longest axis=%.1e, shortest axis=%.1e\n', its, reg, det(Tmve), max(eig(Tmve)), min(eig(Tmve)));
        end

        if its==max_its
            fprintf('Reached the maximum number of iterations, rounding may not be ideal.\n');
        end
    end

    if min(P.b)<=0
        [~,x0] = mve_presolve_cobra(P.A,P.b,150,1e-6);
        [P,N_total,p_shift,T] = shiftPolytope(P,N_total,p_shift,T,eye(dim),x);
        fprintf('Shifting so the origin is inside the polytope...rounding may not be ideal.\n');
    else
        fprintf('Maximum volume ellipsoid found, and the origin is inside the transformed polytope.\n');
    end
    P.p = zeros(dim,1);
else
    P.p = linprog(zeros(size(P.A, 2), 1), P.A, P.b);
    T = eye(size(N_total,2));
end

P.N = N_total;
P.p_shift = p_shift;
P.T = T;
end

function [P,N,p,T] = shiftPolytope(P, N, p, T, trans, shift)
    %shift the polytope by a point and apply a transformation, while retaining
    %the information to undo the transformation later (to recover the samples)
    
    %let x denote the original space, y the current space, and z the new space
    %we have
    %
    %   P.A y <= P.b   and x = N*y+p
    %
    %  applying the transformation
    %
    %   trans * z + shift = y
    %
    % yields the system
    %
    %  x = (N * trans) * z + N*shift + p
    %  (P.A * trans) * z <= P.b - P.A * shift
    p = p + N*shift;
    N = N * trans;
    T = T * trans;
    P.b = P.b - P.A*shift;
    P.A = P.A*trans;
end