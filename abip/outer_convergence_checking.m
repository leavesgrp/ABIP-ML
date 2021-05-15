function [residual, solved, infeasible, unbounded, ttime] = outer_convergence_checking(data, work, barrier, setting, residual, vars, i, j);
    m=vars.m;
    n=vars.n;
    u=vars.u;
    v=vars.v;
    k=vars.k;

    tau = abs(u(end));
    kap = abs(v(end)) / (work.sc_b * work.sc_c * work.scale);
    y   = u(1:m) / tau;
    x   = u(m+1:m+n) / tau;
    s   = v(m+1:m+n) / tau;
    residual.err_pri  = norm(work.D.*(data.A_times(x) - data.b)) / (1 + work.nm_b) / (work.sc_b * work.scale); 
    residual.err_dual = norm(work.E.*(data.AT_times(y) + s - data.c)) / (1 + work.nm_c) / (work.sc_c * work.scale); 
    pobj     = data.c' * x / (work.sc_c * work.sc_b * work.scale);
    dobj     = data.b' * y / (work.sc_c * work.sc_b * work.scale);
    residual.gap      = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
    
    if (data.c'*u(m+1:m+n) < 0)
%         unb_res = norm(E.*data.c) * norm(D.*(data.A * u(m+1:m+n))) / (-data.c'*u(m+1:m+n)) / scale;
        unb_res = norm(work.E.*data.c) * norm(work.D.*(data.A_times(u(m+1:m+n)))) / (-data.c'*u(m+1:m+n)) / work.scale;
    else
        unb_res = inf;
    end
        
    if (-data.b'*u(1:m) < 0)
%         inf_res = norm(D.*data.b) * norm(E.*(data.A' * u(1:m) + v(m+1:m+n))) / (data.b'*u(1:m)) / scale;
        inf_res = norm(work.D.*data.b) * norm(work.E.*(data.AT_times(u(1:m)) + v(m+1:m+n))) / (data.b'*u(1:m)) / work.scale;
    else
        inf_res = inf;
    end
    
    barrier.ratio = barrier.mu/setting.eps;
    residual.error_ratio = max(residual.gap,max(residual.err_pri,residual.err_dual))/setting.eps;
    solved = residual.error_ratio < 1;
    infeasible = inf_res < setting.eps;
    unbounded = unb_res < setting.eps;
        
    ttime = toc;
    
    fprintf('i: %5d, mu: %3.2e, k: %5d presi: %3.7e dresi: %3.7e, dgap: %3.7e, time: %3.2e \n', ...
                i, barrier.mu, k, residual.err_pri, residual.err_dual, residual.gap, ttime);

end