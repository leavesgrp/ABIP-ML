classdef normal_config  < model_interface
    properties
        Q
        A
        b
        c
        m
        n
        sp
    end
    
    methods
        function data = normal_config(dataset)
            temp=dataset;
            data.A=temp.A;
            data.b=temp.b;
            data.c=temp.c;
            [data.m, data.n]=size(data.A);
            data.Q=0;
        end
        
        function params = params_reset(obj, params)
            
        end
        
        function w=A_times(obj, rhs) %
            w=obj.A*rhs;
        end

        function w=AT_times(obj, rhs);  %
             w=obj.A'*rhs;
        end

        function w=Q_times(obj, rhs);  %
            w=obj.Q*rhs;
        end
        
        
        


        function y=solve_lin_sys(obj, w, rhs, setting, init, warm_start, k);     %Ã»¸ÄÍê£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
        % % This function solves the linear system Wx=rhs for
        % % W=[I_2n A;
        % %        A'   -I_4n];
            m=obj.m;
            n=obj.n;
            rho_y=setting.rho_y;

            if w.METHOD==2
                if k == -1 
                    tol = 1e-9 * norm(rhs(1:m)); 
                else
                    tol = max(1e-9, 1e-5 * norm(rhs(1:m)) / (k+1)^2);  
                end

                y               = zeros(m+n, 1); 
                
                [y(1:m), itn]   = obj.pcg_abips(rhs(1:m)+obj.A_times(rhs(m+1:m+n)), warm_start(1:m), w.M, ...
                    rho_y, m, tol, setting.extra_verbose);  
                
                y(m+1:m+n)      = -rhs(m+1:m+n) + obj.AT_times(y(1:m));
            else
                y               = (w.L'\(w.d\(w.L\rhs(w.P))));
                y(w.P)          = y;
            end
        end



        function [x, i] = pcg_abips(obj, b, x, M, rho_x, max_iters, tol, verbose)
            A=obj.A;
            % M is inverse preconditioner
            r   = b - (rho_x*x + A*(A'*x));
            z   = M.*r;
            p   = z;
            ip  = r'*z;
            
            for i = 1: max_iters
                Ap      = rho_x*p + A*(A'*p);
                alpha   = ip/(p'*Ap);
                x       = x + alpha*p;
                r       = r - alpha*Ap;
                resid   = norm(r); 
                if resid < tol
                    if verbose
                        fprintf('CG took %i iterations to converge, resisdual %4f <= tolerance %4f\n', i, resid, tol)
                    end
                    return;
                end
                z       = M.*r;
                ipold   = ip;
                ip      = z'*r;
                beta    = ip/ipold;
                p       = z + beta*p;
            end

            if verbose
                fprintf('CG did not converge within %i iterations, resisdual %4f > tolerance %4f\n', max_iters, resid, tol)
            end
        end


        
        
        function [obj, w] = normalize_data(obj, w, setting)
        scale=setting.scale;
        A     = obj.A;
        b     = obj.b; 
        c     = obj.c; 
        m    = obj.m;
        n     = obj.n; 

        min_scale   = 1e-3;
        max_scale   = 1e3;
        minRowScale = min_scale * sqrt(n);
        maxRowScale = max_scale * sqrt(n);
        minColScale = min_scale * sqrt(m);
        maxColScale = max_scale * sqrt(m);

        %% E scale
        E = sqrt(sum(A.^2, 1))';
        E(E < minColScale) = 1;
        E(E > maxColScale) = maxColScale;
        A = A*sparse(diag(1./E));

        %% D scale:
        D = sqrt(sum(A.^2, 2));
        D(D < minRowScale) = 1;
        D(D > maxRowScale) = maxRowScale;
        A = sparse(diag(1./D))*A;

        %% b and c scale
        nmrowA = mean(sqrt(sum(A.^2, 2)));  
        nmcolA = mean(sqrt(sum(A.^2, 1)));  

        A = A*scale;

        c = c./E;
        cnorm = norm(c); 
        sc_c = nmrowA/max(cnorm, min_scale);
        c = c * sc_c * scale;

        b = b./D;
        bnorm = norm(b);  
        sc_b = nmcolA/ max(bnorm, min_scale);
        b = b * sc_b * scale;



        %% normalized (A,b,c) record
        obj.A = A; 
        obj.b = b;
        obj.c = c;



        w.D = D;
        w.E = E;
        w.sc_b = sc_b;
        w.sc_c = sc_c;

        end
        
        
        function [obj, work]=no_normalize_settings(obj, work, setting);

        end

        
        
        
        function sp=individual_sparsity_ratio(obj);   %
              sp=nnz(obj.A)/obj.m/obj.n;
        end

        
        
        
        
        function [work,obj]=individual_data_work_settings(obj,work,setting,part);
             rho_y=setting.rho_y;

        if part==1        

        elseif part==2
            obj.Q= sparse([zeros(obj.m) obj.A -obj.b; -obj.A' zeros(obj.n) obj.c; obj.b' -obj.c' 0]);
             if work.METHOD==1     %direct solver
                  W                           = sparse([rho_y*speye(obj.m) obj.A; obj.A' -speye(obj.n)]);  
                  [work.L, work.d, work.P]    = ldl(W, 'vector');
             elseif work.METHOD==2           %indirect solver
                  work.M = 1 ./ diag(rho_y*speye(obj.m) + obj.A*obj.A'); % pre-conditioner
             end
        end

        end
        
        
    end
        
   
  
    
    
end