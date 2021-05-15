function [x,y,s,info] = abip(data,params)


import model.*
%% Initial instance
if isfield(params,'Problem')
    data=eval([params.Problem,'_config(data)']);
       fprintf('Problem type: %s\n',params.Problem);
    params = data.params_reset(params);
else
    data=normal_config(data);
    fprintf('Problem type: normal\n');
end




%% default settings
max_outiters = 100; 
max_iters    = 1000000;             
eps          = 1e-3;                   
alpha        = 1.8; 
mu           = 1.0; 
normalize    = 1;                
scale        = 1;  
rho_y        = 1e-3;
sigma        = 0.3; 
adaptive     = 1;                 
eps_cor      = 0.2;
eps_pen      = 0.1; 
% added
SPARSITY_RATIO=0.2;
ADAPTIVE_LOOKBACK=20;
METHOD=1;

% conjugate gradient (CG) settings:
use_indirect    = false;   % use conjugate gradient rather than direct method
extra_verbose   = false;  % CG prints summary

%% constants
undet_tol    = 1e-18;             % tol for undetermined solution (tau = kappa = 0)


%% setting 

if nargin==2
    if isfield(params,'max_ipm_iters');   max_outiters = params.max_ipm_iters;   end
    if isfield(params,'max_admm_iters');  max_iters    = params.max_admm_iters;      end
    if isfield(params,'eps');             eps          = params.eps;            end
    if isfield(params,'alpha');           alpha        = params.alpha;          end
    if isfield(params,'sigma');           sigma        = params.sigma;          end
    if isfield(params,'normalize');       normalize    = params.normalize;      end
    if isfield(params,'scale');           scale        = params.scale;          end
    if isfield(params,'rho_y');           rho_y        = params.rho_y;          end
    if isfield(params,'adaptive');        adaptive     = params.adaptive;       end
    if isfield(params,'eps_cor');         eps_cor      = params.eps_cor;        end  
    if isfield(params,'eps_pen');         eps_pen      = params.eps_pen;        end  
    % added
    if isfield(params,'SPARSITY_RATIO');    SPARSITY_RATIO=params.SPARSITY_RATIO;    end
    if isfield(params,'ADAPTIVE_LOOKBACK');    ADAPTIVE_LOOKBACK=params.ADAPTIVE_LOOKBACK;    end
    if isfield(params,'METHOD');      METHOD=params.METHOD;      end
end

setting=struct();
setting.max_outiters = max_outiters; 
setting.max_iters      = max_iters;             
setting.eps           = eps;                   
setting.alpha        = alpha; 
setting.normalize    = normalize;
setting.scale            =scale;
setting.rho_y        = rho_y;
setting.adaptive     = adaptive;                 
setting.eps_cor      = eps_cor;
setting.eps_pen      = eps_pen; 
setting.SPARSITY_RATIO = SPARSITY_RATIO;
setting.ADAPTIVE_LOOKBACK = ADAPTIVE_LOOKBACK;
setting.METHOD=METHOD;
setting.use_indirect    = false;   % use conjugate gradient rather than direct method
setting.extra_verbose   = false;  % CG prints summary
setting.undet_tol    =    undet_tol;



%% initialization of work, barrier, residual
work=struct();
work.METHOD=METHOD;
barrier=struct();
residual=struct();
[data, work, barrier, residual]=update_work(data, work, setting);



% dimension
n = length(data.c); 
m = length(data.b); 
l = n+m+1;
u = zeros(l, 1);
v = zeros(l, 1);




%% Initialization  (cold start)
u(m+1:l) = ones(n+1,1)*sqrt(barrier.mu/barrier.beta); 
v(m+1:l) = ones(n+1,1)*sqrt(barrier.mu/barrier.beta);
k        = 0; 

vars=struct();
vars.u=u;
vars.v=v;
vars.m=m;
vars.n=n;
vars.l=l;
vars.k=k;


barrier.ratio=barrier.mu/setting.eps;

tic;
for i=0:setting.max_outiters-1   
    % inner_stopper
    if (min(data.sp,setting.SPARSITY_RATIO)>0.5)
        inner_stopper=barrier.mu^(-0.35);
    elseif (min(data.sp,setting.SPARSITY_RATIO)>0.2)
        inner_stopper=barrier.mu^(-1);
    else
        inner_stopper=setting.max_iters;
    end
        
%     inner_stopper=barrier.mu^(-0.5);
    for j=0:inner_stopper

        %% linera_barrier_dual_update
        vars=linear_barrier_dual_update(data, work, setting, barrier, vars);
        
        %% convergence checking: (inner calc_residuals)
        [residual, vars, stat]=inner_convergence_checking(data, work, barrier, setting, residual, vars, j, 1);
        if stat
            break;
        end
        [residual, vars, stat]=inner_convergence_checking(data, work, barrier, setting, residual, vars, j, 2);
        if stat
            break;
        end  
        
    end
    
    %added
    if (vars.k>setting.max_iters*0.8)
        barrier.final_check=1;
    end
    
    %% convergence checking
    [residual, solved, infeasible, unbounded, ttime] = outer_convergence_checking(data, work, barrier, setting, residual, vars, i, j);
    if (solved || infeasible || unbounded)
        break;
    end
    
    if (vars.k+1>max_iters) 
        break;
    end
    
    
    %% update_barrier    
    barrier=update_barrier( residual, data.sp, setting, barrier);
  
    
    %% reinitialization
    vars=reinitialize_vars(barrier.sigma,vars,0);  
    if adaptive
        vars=reinitialize_vars(barrier.sigma,vars,1);
        barrier.beta=1;
        barrier.beta = BBspectral(data, work, barrier, setting, vars);
        vars=reinitialize_vars(barrier.sigma,vars,2);
    end
end



if (vars.k+1 > setting.max_iters); vars.k=vars.k+1; end
if (i+1 == setting.max_outiters); i=i+1; end
ttime = toc;

%% Certificate of infeasibility
[x,y,s,info]=Certificate_of_infeasibility(data, work, setting, residual, vars, ttime, i);

end