function [x,y,s,info]=Certificate_of_infeasibility(data, work, setting, residual, vars, ttime, i);
u=vars.u;
v=vars.v;
ut=vars.ut;
k=vars.k;
m=vars.m;
n=vars.n;
l=m+n+1;



undet_tol=setting.undet_tol;
normalize=setting.normalize;

%% Certificate of infeasibility
tau = abs(u(end));
kap = abs(v(end)) / (work.sc_b * work.sc_c * work.scale);

y = u(1:m) / tau;
x = u(m+1:m+n) / tau;
s = v(m+1:m+n) / tau;

if (tau > undet_tol)
    status = 'solved';
else
    y = nan(m,1);
    x = nan(n,1);
    s = nan(n,1);
    
    y_h = u(1:m);
    x_h = u(m+1:m+n);
    s_h = v(m+1:m+n);
    if norm((u+ut)/2)<=2*(undet_tol*sqrt(l))
        status = 'undetermined';
    elseif data.c'*x_h < data.b'*y_h
        status = 'infeasible';
        y = y_h * work.scale * work.sc_b * work.sc_c /(data.b'*y_h);
        s = s_h * work.scale * work.sc_b * work.sc_c /(data.b'*y_h);
        x = -x_h * work.scale * work.sc_b * work.sc_c /(data.c'*x_h);
    else
        status = 'unbounded';
        y = y_h * work.scale * work.sc_b * work.sc_c /(data.b'*y_h);
        s = s_h * work.scale * work.sc_b * work.sc_c /(data.b'*y_h);
    end
end

info.status    = status;
info.outiter  = i; 
info.iter = k;

info.resPri    = residual.err_pri;
info.resDual   = residual.err_dual;
info.relGap    = residual.gap;
info.time      = ttime; 

if (normalize)
    x = x ./ (work.E * work.sc_b);
    y = y ./ (work.D * work.sc_c);   
    s = s .* (work.E / (work.sc_c * work.scale));
end

info.pobj    = data.c'*x/work.sc_c/setting.scale;  




end