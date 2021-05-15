function [residual, vars, status]=inner_convergence_checking(data, work, barrier, setting, residual, vars, j, part)
    u=vars.u;
    v=vars.v;
    k=vars.k;
    if part==1
        err_inner = norm(data.Q_times(u)-v)/(1+norm([u;v]));
        tol = barrier.gamma*barrier.mu;
        k = k+1;
        vars.k=k;
        status=err_inner < tol;
    elseif part==2
        m=vars.m;
        n=vars.n;
        status=false;
        if (barrier.final_check && mod(j+1,1)==0)
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
            
            residual.error_ratio = max(residual.gap,max(residual.err_pri,residual.err_dual))/setting.eps;
            solved = residual.error_ratio < 1;

            status=(solved || k+1>setting.max_iters);         
        end      
    end
end