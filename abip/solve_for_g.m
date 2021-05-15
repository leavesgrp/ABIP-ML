function [h, g, gTh] = solve_for_g(data, work, setting)
h           = [-data.b; data.c];
g=data.solve_lin_sys(work, h, setting, true, zeros(data.m,1), 0);
g(data.m+1:end)  = -g(data.m+1:end);
gTh         = g'*h;
end