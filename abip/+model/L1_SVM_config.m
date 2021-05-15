classdef L1_SVM_config < model_interface
    properties
        Z
        d
        y
        b
        c
        m
        n
        m_X
        n_X
        sp
        scalar
    end
    
    methods
        function data = L1_SVM_config(dataset)
            temp=dataset;
           %% pre normalize
            X=temp.X;
            F = max(abs(X))';
            X(:,(F==0))=[];
            F(F==0)=[];
            [m,n]=size(X);
            F=max(1e-4,F);
            F= spdiags(1./F,0:0,n,n);
            X=X*F;            
            
            [data.m_X, data.n_X]=size(X);
            data.y=temp.y;
            data.Z=spdiags(data.y,0:0,data.m_X,data.m_X)*X;
            data.sp=nnz(data.Z)/m/n;
            
            data.m=data.m_X;
            data.n=2*data.m_X+2*data.n_X+2;
            data.d=speye(data.m_X);
            
            if isfield(temp,'scalar')
                data.scalar=temp.scalar;
            else
                data.scalar=1;
            end            
        end
        
        function params = params_reset(obj, params)
            if isfield(params,'METHOD') == 0;   params.METHOD = 1;   end
            if isfield(params,'normalize') == 0;   params.normalize = 1;   end
            if isfield(params,'alpha') == 0;   params.alpha = 1.7;   end
        end
        
        function w=A_times(obj, rhs)
            m_X=obj.m_X;
            n_X=obj.n_X;
            w=obj.d*(rhs(1:m_X)-rhs(m_X+1:2*m_X))+obj.Z*(rhs(2*m_X+1:2*m_X+n_X)-rhs(2*m_X+n_X+1:2*m_X+2*n_X))+obj.y*(rhs(2*m_X+2*n_X+1)-rhs(2*m_X+2*n_X+2));
        end

        function w=AT_times(obj, rhs);  %n=n_X
            m_X=obj.m_X;
            n_X=obj.n_X;
            x=obj.d'*rhs;
            y=obj.Z'*rhs;
            z=obj.y'*rhs;
            w=[x;-x;y;-y;z;-z];
        end

        function w=Q_times(obj, rhs);  %n=n_X
            m=obj.m;
            n=obj.n;
            x=obj.A_times(rhs(m+1:m+n))-obj.b*rhs(m+n+1);
            y=-obj.AT_times(rhs(1:m))+obj.c*rhs(m+n+1);
            z=obj.b'*rhs(1:m)-obj.c'*rhs(m+1:m+n);
            w=[x;y;z];
        end
        
        
        function w=solve_lin_sys(obj, work, rhs, setting, init, warm_start, k);
        % % This function solves the linear system Wx=rhs for
        % % W=[I_2n A;
        % %        A'   -I_4n];
        rho_y=setting.rho_y;

        m=obj.m;
        n=obj.n;

        w               = zeros(m+n, 1); 
        w(1:m)=obj.solve_nolmal_equation(work, rhs(1:m)+obj.A_times(rhs(m+1:m+n)), rho_y, setting.use_indirect, init, warm_start, k);
        w(m+1:m+n)      = -rhs(m+1:m+n) + obj.AT_times(w(1:m)) ;

        end

        
        
        
        
        
        
        function w=solve_nolmal_equation(obj, work, rhs, rho_y, use_indirect, init, warm_start, k);
        % this function is to solve 
        % (rho_y*speye(m)+data.A*data.A')\rhs
        if work.METHOD==1
            y=obj.solve_D_plus_yyT(work, rhs);
            z=obj.Z'*y;
            if use_indirect
                fprintf('error!!!!!!!');
            else
                z=work.L\(work.L'\z);    %use direct
            end
            w=0.5*(   y - solve_D_plus_yyT(obj, work,      obj.Z*z       )    );
        elseif work.METHOD==2
            if use_indirect
                fprintf('error!!!!!!!');
            else
                b=rhs*0.5;
                p=work.L\b;
                r=work.v'*p/work.c;
            w=work.L'\(p-work.v*r);
            end
        elseif work.METHOD==3
            w=obj.pcg_abip(work, rhs, rho_y, work.METHOD, init, warm_start, k);
        end
        end
        

        
        
        
        function w=solve_D_plus_yyT(obj, work, rhs);
            w=work.D_inv*rhs-work.c*(obj.y'*(work.D_inv*rhs));
        end


        
        

        function x = pcg_abip(obj, work, rhs, rho_x, METHOD, init, warm_start, k);

            [m,n]=size(obj.Z);
            b=rhs;
            x=warm_start;
            M=work.M;
%             A=obj.A;
            verbose=0;
            % M is inverse preconditioner
            if init 
                tol = 1e-9 * norm(rhs(1:m)); 
            else
                tol = max(1e-9, 1e-5 * norm(rhs(1:m)) / (k+1)^2);  
            end

%             r   = b - (rho_x*x + A*(A'*x));
            r   = b - (rho_x*x + obj.A_times(obj.AT_times(x)));
            z   = M.*r;
            p   = z;
            ip  = r'*z;
            
            max_iters=m;
            for i = 1: max_iters
%                 Ap      = rho_x*p + A*(A'*p);
                Ap      = rho_x*p + obj.A_times(obj.AT_times(p));
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
        b     = obj.b; 
        c     = obj.c; 
        m=obj.m;
        n=obj.n;
        m_X=obj.m_X;
        n_X=obj.n_X;

        min_scale   = 1e-3;
        max_scale   = 1e3;
        minRowScale = min_scale * sqrt(n);
        maxRowScale = max_scale * sqrt(n);
        minColScale = min_scale * sqrt(m);
        maxColScale = max_scale * sqrt(m);

        %% E scale
        nor_scalar=1;
        sp_Z=nnz(obj.Z)/obj.m_X/obj.n_X; %
        if sp_Z<0.65 
            nor_scalar=0.15;
            if sp_Z<0.01
                nor_scalar=0.9;
            end
        end
        
        x1=nor_scalar*sqrt(sum(obj.Z.^2, 1))';
        x2=nor_scalar*sqrt(sum(obj.y.^2, 1))';
        E = [ones(m_X,1); ones(m_X,1); x1; x1; x2; x2];
        E(E < minColScale) = 1;
        E(E > maxColScale) = maxColScale;
        obj.Z=obj.Z*spdiags(1./x1, 0:0, n_X,n_X);
        obj.y=obj.y/x2;

        %% D scale:
        D=2*(sum(obj.Z.^2, 2)+obj.y.^2+1);
        D=sqrt(D);
        D(D < minRowScale) = 1;
        D(D > maxRowScale) = maxRowScale;
        obj.d=sparse(spdiags(1./D, 0:0, m,m));
        obj.Z=obj.d*obj.Z;
        obj.y=obj.d*obj.y;

        %% b and c scale
        nmrowA=mean(sqrt(2*(sum(obj.Z.^2, 2)+obj.y.^2+sum(obj.d.^2, 2))));
        x1=sqrt(sum(obj.d.^2, 1))';
        x2=sqrt(sum(obj.Z.^2, 1))';
        x3=sqrt(sum(obj.y.^2, 1))';
        nmcolA = mean([x1;x1;x2;x2;x3;x3]);  

        obj.d=obj.d*scale;                 %change d, Z, y
        obj.Z=obj.Z*scale;
        obj.y=obj.y*scale;

        c = c./E;
        cnorm = norm(c); 
        sc_c = nmrowA/max(cnorm, min_scale);
        c = c * sc_c * scale;

        b = b./D;
        bnorm = norm(b); 
        sc_b = nmcolA/ max(bnorm, min_scale);
        b = b * sc_b * scale;



        %% normalized (A,b,c) record
        obj.b = b;
        obj.c = c;

        w.D = D;
        w.E = E;
        w.sc_b = sc_b;
        w.sc_c = sc_c;

        end
        
        
        
        function [obj, work]=no_normalize_settings(obj, work, setting);

        end


        function sp=individual_sparsity_ratio(obj);
              sp=2*(2*obj.m_X+nnz(obj.Z))/obj.m/obj.n;
        end




        function [work,obj]=individual_data_work_settings(obj,work,setting,part);
             rho_y=setting.rho_y;
             m_X=obj.m_X;
             n_X=obj.n_X;
             m=obj.m;
             n=obj.n;

        if part==1
            work.lambda=sqrt(m_X*log(n_X));
            work.lambda=work.lambda*obj.scalar;
            obj.b=ones(m_X,1);
            obj.c=[ones(m_X,1);zeros(m_X,1);work.lambda*ones(n_X,1);work.lambda*ones(n_X,1);0;0];
        elseif part==2
              sp_Z=nnz(obj.Z)/m_X/n_X;                     % use sparse format or dense format
              if sp_Z>0.8
                  obj.Z=full(obj.Z);
              end

              tic;
              if work.METHOD==1
                  D=obj.d*obj.d+0.5*rho_y*speye(m_X);                   
                  work.D_inv=inv(D);
                  work.L=full(chol(      speye(n_X)  +  obj.Z'*work.D_inv*obj.Z  -  (obj.Z'*(work.D_inv*obj.y))*(1+obj.y'*work.D_inv*obj.y)^(-1)*(obj.y'*work.D_inv*obj.Z)     ));
                  work.c=(1+obj.y'*work.D_inv*obj.y)^(-1)*(work.D_inv*obj.y);
                  % time and condition
                  sol_time = toc; 
                  L=work.L;
                  b=whos('L');
                  sparsity_L=nnz(work.L)/(n_X^2);
                  fprintf('Matrix decomposition:  %3.2e \n', sol_time);
                  fprintf('memory of L:  %3.2e MB,  sparsity of L: %3.2e \n', b.bytes/1024/1024, sparsity_L);
              elseif work.METHOD==2
                  work.L=(sparse(chol(            obj.d*obj.d+0.5*rho_y*speye(m_X)       +obj.Z*obj.Z'               )))';   %LL'=W
                  work.v=work.L\obj.y;
                  work.c=work.v'*work.v+1;
                  % time and condition
                  sol_time = toc; 
                  L=work.L;
                  b=whos('L');
                  sparsity_L=nnz(work.L)/(n_X^2);
                  fprintf('Matrix decomposition:  %3.2e \n', sol_time);
                  fprintf('memory of L:  %3.2e MB,  sparsity of L: %3.2e \n', b.bytes/1024/1024, sparsity_L);
              elseif work.METHOD==3
                  work.M = 1 ./ diag(rho_y*speye(m) + obj.A*obj.A'); % pre-conditioner
              end

        end



        end

        
 
    end
 
    
    
end


