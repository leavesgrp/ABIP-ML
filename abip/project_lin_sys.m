function ut = project_lin_sys(data, work, setting, k, u, v)    
    rho_y=setting.rho_y;

    m=data.m;
    n=data.n;
    
    h=work.h;
    g=work.g;
    gTh=work.gTh;
    
    ut                  = u+v;
    ut(1:m)             = rho_y*ut(1:m);
    ut(1:m+n)           = ut(1:m+n) - ut(end)*h;
    ut(1:m+n)           = ut(1:m+n) - h*((g'*ut(1:m+n))/(gTh+1));
    warm_start          = u(1:n+m);
    ut(m+1:end-1)       = -ut(m+1:end-1);
    ut(1:m+n)    = data.solve_lin_sys(work, ut(1:m+n), setting, false, warm_start(1:m), k);
    ut(end)             = (ut(end) + h'*ut(1:m+n));   
end