function [x4, v4, step, done_ReverseCheck] = generalized_leapfrog(x0, v0, h, ham, opts)
    % Step 1
    x1 = x0;
    v1 = v0 - (h/2) * ham.DU(x0);
    done = 0;
    done_ReverseCheck = 1;
    
    % Step 2
    x2 = x1; v2 = v1;
    nu = zeros(size(x1,1),size(ham.A,1),1);
    [dKdv, dKdx, nu] = ham.approxDK(x1, v2, nu);
    dKdv_x2 = dKdv;
    for step = 1:opts.maxODEStep
        x2_old = x2;
        v2 = v1 - h * dKdx/2;
        x2 = x1 + h * (dKdv+dKdv_x2)/2;
        dist = ham.x_norm(x2_old, x2-x2_old)./ h;
        if (max(dist,[],'all') < opts.implicitTol)
            done = 1;
            break;
        elseif any(dist > 1e16, 'all')
            break;
        end
        [dKdv, dKdx, nu] = ham.approxDK(x1, v2, nu);
        dKdv_x2 = ham.approxDKv(x2, v2, nu);
    end
    [~, dKdx, ~] = ham.approxDK(x2, v2, nu);    
    
    if done == 0
        x4 = NaN; v4 = NaN;
        return
    else
        v2 = v2 - h * dKdx/2;
    end
    
    done_back = 0;
    x1_back = x2;
    v1_back = -v2;
    v2_back = -v2;
    x2_back = x2;
    nu_back = zeros(size(x1_back,1),size(ham.A,1),1);
    [dKdv_back, dKdx_back, nu_back] = ham.approxDK(x1_back, v2_back, nu_back);
    dKdv_x2_back = dKdv_back;
    for step = 1:opts.maxODEStep
        x2_old_back = x2_back;
        v2_back = v1_back - h * dKdx_back/2;
        x2_back = x1_back + h * (dKdv_back+dKdv_x2_back)/2;
        dist = ham.x_norm(x2_old_back, x2_back-x2_old_back)./ h;
        if (max(dist,[],'all') < opts.implicitTol)
            done_back = 1;
            break;
        elseif any(dist > 1e16, 'all')
            break;
        end
        [dKdv_back, dKdx_back, nu_back] = ham.approxDK(x1_back, v2_back, nu_back);
        dKdv_x2_back = ham.approxDKv(x2_back, v2_back, nu_back);
    end
    [~, dKdx_back, ~] = ham.approxDK(x2_back, v2_back, nu_back);    
    
    if done_back == 0
        x4 = NaN; v4 = NaN;
        return
    else
        v2_back = v2_back - h * dKdx_back/2;
    end
    
    v2_back = - v2_back;
    dist = ham.x_norm(x1, x2_back-x1) + ham.x_norm(x2_back, x2_back-x1) + ham.v_norm(v1, v2_back-v1) + ham.v_norm(v2_back, v2_back-v1);
    done_ReverseCheck = dist < opts.implicitTol_RC;
    %x4 = x0; v4=v0;
    
    % Step 3
    x3 = x2;
    v3 = v2 - (h/2) * ham.DU(x3);
    
    % Step 4 (Project to Ax = b)
    ham.prepare(x3);
    v4 = v3;
    x4 = ham.project(x3);
end
