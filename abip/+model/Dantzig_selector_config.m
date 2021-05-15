classdef Dantzig_selector_config < model_interface
    properties
        X
        y
        b
        c
        m
        n
        m_X
        n_X
        sp
        scalar
        normalize
        ld
        rd
        md
        d
        dI
        dI_inv
        ld_invl
    end
    
    methods
        function data = Dantzig_selector_config(dataset)
            temp=dataset;
           %% pre normalize
            X=temp.X;
            
           
            
            F=sqrt(sum(X.^2,1))';
            X(:,(F==0))=[];
            F(F==0)=[];
            [m,n]=size(X);
            F=max(1e-4,F);
            F= spdiags(1./F,0:0,n,n);
            X=X*F; 
            
            sp_X=nnz(X)/m/n;
            if sp_X>0.5
                X=full(X);
            end
            
            data.X=X;
            data.y=temp.y;
            data.b=0;
            data.c=0;
            data.m=0;
            data.n=0;
            data.m_X=0;
            data.n_X=0;
            if isfield(temp,'scalar')
                data.scalar=temp.scalar;
            else
                data.scalar=1;
            end
        end
        
        function params = params_reset(obj, params)
            if isfield(params,'METHOD') == 0;   params.METHOD = 1;   end
            if isfield(params,'normalize') == 0;   params.normalize = 0;   end
            if isfield(params,'alpha') == 0;   params.alpha = 1.8;   end
        end
        
        function w=A_times(obj, rhs)
        if obj.normalize==1 
            m=obj.m;
            n=obj.n;
            n_X=obj.n_X;           
            s=obj.ld*(   obj.X'*(       obj.X*(        obj.rd*(      rhs(m+1:m+n_X)-rhs(m+n_X+1:n)     )          )        )       );
            w=[obj.d*rhs(1:n_X)+s; obj.d*rhs(n_X+1:m)-s];
        else
            n=obj.n_X;
            y=obj.X'*(obj.X*(rhs(2*n+1:3*n)-rhs(3*n+1:4*n)));
            w=rhs(1:2*n)+[y;-y]; 
        end
        end

        function w=AT_times(obj, rhs);
        if obj.normalize==1
            s=obj.rd*(  obj.X'*(  obj.X*(   obj.ld*(    rhs(1:obj.n_X)-rhs(obj.n_X+1:obj.m)       )     )     )    );
            w=[obj.d*rhs(1:obj.n_X) ; obj.d*rhs(obj.n_X+1:obj.m) ; s ; -s  ];
        else
            n=obj.n_X;
            y=obj.X'*(obj.X*(rhs(1:n)-rhs(n+1:2*n)));
            w=[rhs(1:2*n);y;-y];            
        end
        end

        
        
        function w=Q_times(obj, rhs);
            m=obj.m;
            n=obj.n;
            x=obj.A_times(rhs(m+1:m+n))-obj.b*rhs(m+n+1);
            y=-obj.AT_times(rhs(1:m))+obj.c*rhs(m+n+1);
            z=obj.b'*rhs(1:m)-obj.c'*rhs(m+1:m+n);
            w=[x;y;z];
        end
        
        
        function w=XXT_times(obj, work, rhs)
            %according to sparsity
            if work.use_XXT==1;
                w=work.XXT*rhs;
            else
                w=obj.X*(   obj.X'*rhs   );
            end
        end        
        
        
        
        
        function w=XmdXT_times(obj, work, rhs)
            if work.use_XmdXT==1;
                w=work.XmdXT*rhs;
            else
                w=obj.X*(  obj.md*( obj.X'*rhs )  );                
            end
        end
        
        function w=ldXTXmdX_times(obj, work, rhs)
            w=obj.ld*(     obj.X'*(     obj.XmdXT_times(work, rhs)          )       );
        end
       
        
        
        function w=solve_lin_sys(obj, work, rhs, setting, init, warm_start, k);
            % % This function solves the linear system Wx=rhs for
            % % W=[I_2n A;
            % %        A'   -I_4n];
            rho_y=setting.rho_y;
            m=obj.m;
            n=obj.n;
            
            
            tol=1;
            if work.METHOD==3
                if k == -1 
                    tol = 1e-10 * norm(rhs(1:m)); 
                else
                    tol = max(1e-10, 1e-7 * norm(rhs(1:m)) / (k+1)^2);  
                end   
            end
                 
            
            w               = zeros(m+n, 1); 
            w(1:m)=obj.solve_nolmal_equation(work, rhs(1:m)+obj.A_times(rhs(m+1:m+n)), rho_y, init, warm_start, tol, k);
            w(m+1:m+n)      = -rhs(m+1:m+n) + obj.AT_times(w(1:m)) ;
        end
          

        function w=I_plus_cXXTXXT(obj, work, rho_y, rhs)
            w=rhs+4/(rho_y+1)*obj.XXT_times(work, obj.XXT_times(work, rhs)  );
        end 
        
        
        function w=solve_nolmal_equation(obj, work, rhs, rho_y, init, warm_start, tol, k);
            % this function is to solve 
            % (rho_y*speye(m)+data.A*data.A')\rhs
        if obj.normalize==1
            if work.METHOD==1
                m_X=obj.m_X;
                n_X=obj.n_X;
                m=obj.m;
                n=obj.n;

                b=0.5*rhs;
                z=[obj.dI_inv*b(1:n_X); obj.dI_inv*b(n_X+1:m)];
                s=      obj.X*(      obj.ld*(   z(1:n_X)-z(n_X+1:m)   )             )  ;
                s=obj.ldXTXmdX_times(  work,  work.U\(work.L\s)    ); 
                w=z-[obj.dI_inv*s;-obj.dI_inv*s];
            elseif work.METHOD==3
      
                

            end
        else
            if work.METHOD==1
                b=0.5*rhs;
                kappa=2/(1+rho_y);
                z=kappa*b;
                s=      obj.X*(z(1:obj.n_X)-z(obj.n_X+1:2*obj.n_X))               ;
                s=obj.X'*(  obj.XXT_times(  work,  work.L\(work.L'\s)    )          ); 
                w=z-kappa*[s;-s];
            elseif work.METHOD==3
                
                m_X=obj.m_X;
                n_X=obj.n_X;
                m=obj.m;
                n=obj.n;
                              
                b=0.5*warm_start;
                kappa=2/(1+rho_y);
                z=kappa*b;
                warm_start=      obj.X*(z(1:obj.n_X)-z(obj.n_X+1:2*obj.n_X))               ;
                
                % b
                b=0.5*rhs;
                kappa=2/(1+rho_y);
                z=kappa*b;
                s=      obj.X*(z(1:obj.n_X)-z(obj.n_X+1:2*obj.n_X))   ;
                
                % pcg
                max_iters=m_X;
                verbose=0;
                pcg_sol = pcg_abips(obj, work, s, warm_start, work.M, rho_y, max_iters, tol, verbose);
                
                %pcg_solution
                s=obj.X'*(  obj.XXT_times(  work,  pcg_sol    )          ); 
                w=z-kappa*[s;-s];   

            end
        end
        end
            

        
        
        
        function x = pcg_abips(obj, work, b, x, M, rho_x, max_iters, tol, verbose)
            % A=obj.A;
            % M is inverse preconditioner
            % r   = b - (rho_x*x + A*(A'*x));
            r   = b-obj.I_plus_cXXTXXT(work, rho_x, x);
            z   = M.*r;
            p   = z;
            ip  = r'*z;
            
            for i = 1: max_iters
%                 Ap      = rho_x*p + A*(A'*p);
                Ap       = obj.I_plus_cXXTXXT(work, rho_x, p);
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
        
        
        
        
        
        
        
        
        
        

        
        
        function [obj, work] = normalize_data(obj, work, setting)
            scale=setting.scale;
%             A     = obj.A;
            b     = obj.b; 
            c     = obj.c; 
%             [m,n] = size(A);
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
            
            % if the size is not very big, we can directly calculate XXT
            XTX=obj.X'*obj.X;
            
            %% E scale
            XE=sum(XTX.^2, 1)';
            E=[ ones(m,1); sqrt(2*XE); sqrt(2*XE) ];
            E(E < minColScale) = 1;
            E(E > maxColScale) = maxColScale;
            obj.rd=spdiags( 1./E(m+1:m+n_X)  ,0:0,  n_X,n_X  );

            %% D scale:
            XD=sum(     (    XTX*obj.rd    ).^2, 2);
            D = [ sqrt(1+2*XD)  ;  sqrt(1+2*XD)  ];
            D(D < minRowScale) = 1;
            D(D > maxRowScale) = maxRowScale;
            obj.ld= spdiags( 1./D(1:n_X)  ,0:0,  n_X,n_X  );
            
            obj.d=obj.ld;
            obj.md=obj.rd*obj.rd;
 
            %% b and c scale
%             nmrowA = mean(sqrt(sum(A.^2, 2)));
            XD=sum(  (obj.ld*XTX*obj.rd).^2, 2   );
            nmrowA=mean( sqrt([2*XD+diag(obj.d).^2; 2*XD+diag(obj.d).^2  ]) );
%             nmcolA = mean(sqrt(sum(A.^2, 1)));  
            XE=sum(  (obj.ld*XTX*obj.rd).^2, 1   )';
            nmcolA=mean( sqrt( [diag(obj.d).^2; diag(obj.d).^2 ; 2*XE; 2*XE  ]  ) );

%             A = A*scale;

            c = c./E;
            cnorm = norm(c); 
            % cnorm = sum(abs(c));
            sc_c = nmrowA/max(cnorm, min_scale);
            c = c * sc_c * scale;

            b = b./D;
            bnorm = norm(b); 
            % bnorm = sum(abs(b)); 
            sc_b = nmcolA/ max(bnorm, min_scale);
            b = b * sc_b * scale;



            %% normalized (A,b,c) record
            obj.b = b;
            obj.c = c;

            work.D = D;
            work.E = E;
            work.sc_b = sc_b;
            work.sc_c = sc_c;

        end
        
        
        function [obj, work]=no_normalize_settings(obj, work, setting);
            obj.rd=speye(obj.n_X);
            obj.ld=speye(obj.n_X);
            obj.md=speye(obj.n_X);
            obj.d=speye(obj.n_X);
        end

        function sp=individual_sparsity_ratio(obj);
              sp=0.5;
        end

        
        function [work,obj]=individual_data_work_settings(obj,work,setting,part);
             rho_y=setting.rho_y;
             [m_X,n_X]=size(obj.X); 

        if part==1      
             lambda=sqrt(log(n_X));
             lambda=lambda*obj.scalar;
             obj.b=[obj.X'*obj.y;-obj.X'*obj.y]+lambda*ones(2*n_X,1); 
             obj.c=[zeros(2*n_X,1);ones(2*n_X,1)];
             obj.m_X=m_X;
             obj.n_X=n_X;
             obj.m=2*n_X;
             obj.n=4*n_X;
             obj.normalize=setting.normalize;
        elseif part==2
            obj.dI=0.5*rho_y*speye(n_X)+0.5*obj.d*obj.d;
            obj.dI_inv=spdiags( 1./diag(obj.dI) ,0:0, n_X,n_X   );
            % determine use XXT
            nnz_X=nnz(obj.X);
            if 2*nnz_X>obj.m_X^2
                work.use_XmdXT=1;
                work.use_XXT=1;
            else
                work.use_XmdXT=0;
                work.use_XXT=0;
            end
            
            if obj.normalize==1
            
                 if work.METHOD==1
                      fprintf('start calculating XXT\n');
                      tic;
                      C=(obj.X*obj.md)*obj.X';
                      if work.use_XmdXT==1
                          work.XmdXT=C;
                      end
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                      fprintf('start calculate (XTX)2\n');
                      tic;
                      obj.ld_invl=obj.ld*obj.dI_inv*obj.ld;
                      C=speye(obj.m_X)+2*( obj.X*( obj.ld_invl*(  obj.X'*C   )       )      ); 
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                      fprintf('start cholesky\n');
                      tic;
                      [work.L, work.U]=lu(C);             
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                 elseif work.METHOD==3           %general indirect solver
                     work.M=ones(2*n_X,1);
                 end
            else
             if work.METHOD==1
                  fprintf('start calculating XXT\n');
                  tic;
                  C=full(obj.X*obj.X');
                  if work.use_XXT==1
                      work.XXT=C;
                  end
                  time_fac=toc;
                  fprintf('time: %3.2e\n', time_fac);
                  fprintf('start calculate (XTX)2\n');
                  tic;
                  R=speye(obj.m_X)+4/(1+rho_y)*(C'*C); 
                  time_fac=toc;
                  fprintf('time: %3.2e\n', time_fac);
                  fprintf('start cholesky\n');
                  tic;
                  work.L=chol(R);              %  ...=L'*L
                  time_fac=toc;
                  fprintf('time: %3.2e\n', time_fac);
             elseif work.METHOD==3           %general indirect solver
                 work.M=ones(m_X,1);
                 work.use_XXT=0;
             end                
            end

        end

        end

        
        
    end
        
   
  
    
    
end


