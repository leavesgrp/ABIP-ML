function [vars]=linear_barrier_dual_update(data, work, setting, barrier, vars);
        m=vars.m;
        n=vars.n;
        k=vars.k;
        u=vars.u;
        v=vars.v;
        %% solve linear system
        ut = project_lin_sys(data, work, setting, k, u, v);
        %% solve barrier subproblem
        rel_ut      = setting.alpha*ut+(1-setting.alpha)*u;
        rel_ut(1:m) = ut(1:m);                       
        u           = rel_ut - v;
        temp        = u(m+1:end)/2;
        u(m+1:end)  = temp+sqrt(temp.*temp+barrier.mu/barrier.beta);
        %% dual update:
        v = v + (u - rel_ut);
        
        vars.u=u;
        vars.v=v;
        vars.ut=ut;
        vars.k=k;
        

end