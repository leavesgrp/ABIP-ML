function [data, work, barrier, residual]=update_work(data,work,setting)

      SPARSITY_RATIO = setting.SPARSITY_RATIO;
              use_indirect = setting.use_indirect;
                 normalize = setting.normalize;
                        rho_y = setting.rho_y;
            extra_verbose = setting.extra_verbose;
                         scale = setting.scale;
      %% set data.b,data.c
      [work,data]=data.individual_data_work_settings(work,setting,1);
      
      %% set data.m, data.n
      n = length(data.c); 
      m = length(data.b); 
      
      %% set data.sp
      sp=data.individual_sparsity_ratio();
      data.sp=sp;
      
      %% calculate
      nm_c = norm(data.c); 
      nm_b = norm(data.b);

      %% data normalization 
      if (normalize)
          [data, work] = data.normalize_data(work, setting);
          D = work.D;
          E = work.E;
          sc_b = work.sc_b;
          sc_c = work.sc_c;
      else
          [data, work] = data.no_normalize_settings(work, setting); 
          scale = 1;
          D = ones(data.m,1);
          E = ones(data.n,1);
          sc_b = 1;
          sc_c = 1;
      end
      %% initial work
     [work,data]=data.individual_data_work_settings(work,setting,2);

      %% new sigma&gamma method
      if (max(sp,SPARSITY_RATIO)>0.4 || min(sp,SPARSITY_RATIO)>0.1 && min(sp,SPARSITY_RATIO)<0.2)
          sigma=0.3;
          gamma=2.0;
      elseif (min(sp,SPARSITY_RATIO)>0.2)
          sigma=0.5;
          gamma=3.0;
      else
          sigma=0.8;
          gamma=3.0;
      end
       
       %% calculate h g g'h
      [h, g, gTh] = solve_for_g(data, work, setting);

      final_check = 0;
      double_check=0;
      mu=1;
      beta = 1;

     
      %% normal work settings
      work.nm_b=nm_b;
      work.nm_c=nm_c;
      work.scale=scale;
      work.D=D;
      work.E=E;
      work.sc_b=sc_b;
      work.sc_c=sc_c;
      work.h=h;
      work.g=g;
      work.gTh=gTh;

      %% normal barrier settings
      barrier.sigma=sigma;
      barrier.gamma=gamma;
      barrier.mu=mu;
      barrier.final_check=final_check;
      barrier.double_check=double_check;
      barrier.beta=beta;
      
      
      residual.err_pri = 0; 
      residual.err_dual = 0;
      residual.gap = 0;
      
end