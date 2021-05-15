classdef CLIME_config < model_interface
    properties
        X
        EX
        b
        c
        m
        n
        n_X
        d_X
        problem_col
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
        function data = CLIME_config(dataset)
            temp=dataset;
           %% pre normalize
            X=temp.X;

            F1=mean(X);
            [m,n]=size(X);
            F=sum(X.^2)-2*F1.*sum(X)+m*F1.^2;
            F=F'/sqrt(m);

            X(:,(F==0))=[];
            F(F==0)=[];
            [m,n]=size(X);
            F=max(1e-4,F);
            F= spdiags(1./F,0:0,n,n);
            X=X*F; 
            
            % determine whether X is dense or sparse
            [m,n]=size(X);
            sp_X=nnz(X)/m/n;
            if sp_X>0.6
                X=full(X);
            end
            
            data.X=X;
            [n_X,d_X]=size(data.X);
            data.EX=sum(data.X,1)/n_X;
            data.problem_col=temp.problem_col;
            data.b=0;
            data.c=0;
            data.n_X=n_X;
            data.d_X=d_X;            
            data.m=2*d_X;
            data.n=4*d_X;
            if isfield(temp,'scalar')
                data.scalar=temp.scalar;
            else
                data.scalar=1;
            end
        end
        
        function w=Z_times(obj, rhs)
            w=1/sqrt(obj.n_X-1)*(      obj.X*rhs   -  ones(obj.n_X,1)*(      obj.EX*rhs         )                     );
        end
        
        function w=ZT_times(obj, rhs)
            w=1/sqrt(obj.n_X-1)*(      obj.X'*rhs   -  obj.EX'*(      ones(1, obj.n_X)*rhs         )                     );
        end
        

        
            
        
        function w=A_times(obj, rhs)
        if obj.normalize==1
            m=obj.m;
            n=obj.n;
            d_X=obj.d_X;
            s=obj.ld*(   obj.ZT_times(       obj.Z_times(        obj.rd*(      rhs(m+1:m+d_X)-rhs(m+d_X+1:n)     )          )        )       );
            w=[obj.d*rhs(1:d_X)+s; obj.d*rhs(d_X+1:m)-s];
        else
            n=obj.d_X;
            y=obj.ZT_times(            obj.Z_times(    rhs(2*n+1:3*n)-rhs(3*n+1:4*n)    )            );
            w=rhs(1:2*n)+[y;-y];                
        end
        end

        function w=AT_times(obj, rhs);  %n=n_X
        if obj.normalize==1
            s=obj.rd*(  obj.ZT_times(  obj.Z_times(   obj.ld*(    rhs(1:obj.d_X)-rhs(obj.d_X+1:obj.m)       )     )     )    );
            w=[obj.d*rhs(1:obj.d_X) ; obj.d*rhs(obj.d_X+1:obj.m) ; s ; -s  ];
        else
            n=obj.d_X;
            y=obj.ZT_times(obj.Z_times(rhs(1:n)-rhs(n+1:2*n)));
            w=[rhs(1:2*n);y;-y];            
        end
        end

        
        
        function w=Q_times(obj, rhs);  %n=n_X
            m=obj.m;
            n=obj.n;
            x=obj.A_times(rhs(m+1:m+n))-obj.b*rhs(m+n+1);
            y=-obj.AT_times(rhs(1:m))+obj.c*rhs(m+n+1);
            z=obj.b'*rhs(1:m)-obj.c'*rhs(m+1:m+n);
            w=[x;y;z];
        end
        
        
         function w=ZZT_times(obj, work, rhs)
            if work.use_ZZT==1
                w=work.ZZT*rhs;
            else
                w=obj.Z_times(    obj.ZT_times(  rhs )    );
            end
         end       
        
         
        function w=ZmdZT_times(obj, work, rhs)
            % according to sparsity
            if work.use_ZmdZT==1;
                w=work.ZmdZT*rhs;
            else
                w=obj.Z_times(  obj.md*( obj.ZT_times(rhs) )  );                
            end
        end        

        function w=ldZTZmdZ_times(obj, work, rhs)
            % according to sparsity
            w=obj.ld*(     obj.ZT_times(     obj.ZmdZT_times(work, rhs)          )       );
        end
        
        
        function w=solve_lin_sys(obj, work, rhs, setting, init, warm_start, k);
            % % This function solves the linear system Wx=rhs for
            % % W=[I_2n A;
            % %        A'   -I_4n];
            rho_y=setting.rho_y;

            m=obj.m;
            n=obj.n;
            
            w               = zeros(m+n, 1); 
            w(1:m)=obj.solve_nolmal_equation(work, rhs(1:m)+obj.A_times(rhs(m+1:m+n)), rho_y, init, warm_start, k);
            w(m+1:m+n)      = -rhs(m+1:m+n) + obj.AT_times(w(1:m)) ;
        end
          
        function w=solve_nolmal_equation(obj, work, rhs, rho_y, init, warm_start, k);
            % this function is to solve 
            % (rho_y*speye(m)+data.A*data.A')\rhs            
        if obj.normalize==1
            if work.METHOD==1
                d_X=obj.d_X;
                n_X=obj.n_X;
                m=obj.m;
                n=obj.n;

                b=0.5*rhs;
                z=[obj.dI_inv*b(1:d_X); obj.dI_inv*b(d_X+1:m)];
                s=      obj.Z_times(      obj.ld*(   z(1:d_X)-z(d_X+1:m)   )             )  ;
                s=obj.ldZTZmdZ_times(  work,  work.U\(work.L\s)    ); 
                w=z-[obj.dI_inv*s;-obj.dI_inv*s];
            elseif work.METHOD==3 

            end
        else
            if work.METHOD==1
                b=0.5*rhs;
                kappa=2/(1+rho_y);
                z=kappa*b;
                s=  obj.Z_times(z(1:obj.d_X)-z(obj.d_X+1:2*obj.d_X))               ;
                s=obj.ZT_times(         obj.ZZT_times(   work,  (work.L\(work.L'\s))       )       );
                w=z-kappa*[s;-s];
            elseif work.METHOD==3  %

            end           
        end
        end

        function [obj, work] = normalize_data(obj, work, setting) 
            sc=1;
            obj.b=sc*obj.b;
            obj.c=sc*obj.c;
            obj.normalize=0;
            work.sc_b=sc;
            work.sc_c=sc;
            work.D=ones(obj.m,1);
            work.E=ones(obj.n,1);
            
            obj.rd=speye(obj.d_X);
            obj.ld=speye(obj.d_X);
            obj.md=speye(obj.d_X);
            obj.d=speye(obj.d_X);
            obj.b=sqrt(obj.n_X)*obj.b;
        end        
        
        
        
        
        function [obj, work]=no_normalize_settings(obj, work, setting);  % ÐÂ¼ÓµÄ£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
            obj.rd=speye(obj.d_X);
            obj.ld=speye(obj.d_X);
            obj.md=speye(obj.d_X);
            obj.d=speye(obj.d_X);
            obj.b=sqrt(obj.n_X)*obj.b;
        end

        function sp=individual_sparsity_ratio(obj);
              sp=0.5;
        end

        
        function [work,obj]=individual_data_work_settings(obj,work,setting,part);
             rho_y=setting.rho_y;
             m=obj.m;
             n=obj.n;
             d_X=obj.d_X;
             n_X=obj.n_X;

        if part==1      
            lambda=1;
            lambda=lambda*obj.scalar;
            ei=zeros(obj.d_X,1);
            ei(obj.problem_col)=1;
            obj.b=lambda*ones(obj.m,1)+[ei;-ei];
            obj.c=[zeros(2*obj.d_X,1) ; ones(2*obj.d_X,1)];
            obj.normalize=setting.normalize;
        elseif part==2
            obj.dI=0.5*rho_y*speye(d_X)+0.5*obj.d*obj.d;
            obj.dI_inv=spdiags( 1./diag(obj.dI) ,0:0, d_X,d_X   );
            % determine use ZZT
            nnz_Z=nnz(obj.X);
            if 2*nnz_Z>obj.n_X^2
                work.use_ZmdZT=1;
                work.use_ZZT=1;
            else
                work.use_ZmdZT=0;
                work.use_ZZT=1;
            end
            
            if obj.normalize==1
                 if work.METHOD==1
                      fprintf('start calculating ZZT\n');
                      tic;
                          C=1/(n_X-1)*(obj.X-ones(n_X,1)*obj.EX)*obj.md*(obj.X-ones(n_X,1)*obj.EX)';
                      if work.use_ZmdZT==1   
                          work.ZmdZT=full(C);
                      end
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                      fprintf('start calculating (ZZT)^2\n');
                      tic;
                      obj.ld_invl=obj.ld*obj.dI_inv*obj.ld;
                      C=speye(obj.n_X)+2*( obj.Z_times( obj.ld_invl*(  obj.ZT_times(C)   )       )      ); 
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                      fprintf('start cholesky\n');
                      tic;
                      [work.L, work.U]=lu(C);           
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                 elseif work.METHOD==3           %general indirect solver
                     work.M=ones(2*d_X,1);

                 end
            else
                 if work.METHOD==1
                      fprintf('start calculating ZZT\n');
                      tic;
                      C=ones(obj.n_X,1)*( obj.EX*obj.X' );
                      C=1/(obj.n_X-1)*(  obj.X*obj.X' - C-C' + obj.EX*obj.EX'*ones(obj.n_X)            );
                      C=full(C);
                      if work.use_ZZT==1
                          work.ZZT=C;
                      end
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                      fprintf('start calculating (ZZT)^2\n');
                      tic;
                      R=speye(obj.n_X)+4/(1+rho_y)*(C'*C); 
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                      fprintf('start cholesky\n');
                      tic;
                      work.L=chol(R);              %  ...=L'*L
                      time_fac=toc;
                      fprintf('time: %3.2e\n', time_fac);
                 elseif work.METHOD==3           %general indirect solver
                     work.M=ones(2*d_X,1);

                 end                
            end

        end



        end

        
        
    end
        
   
  
    
    
end


