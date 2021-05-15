function barrier=update_barrier(residual,sp,setting, barrier)
      
      SPARSITY_RATIO=setting.SPARSITY_RATIO;

      error_ratio=residual.error_ratio;

      sigma=barrier.sigma;
      gamma=barrier.gamma;
      mu=barrier.mu;
      final_check=barrier.final_check;
      double_check=barrier.double_check;
      ratio=barrier.ratio;

if (max(sp,SPARSITY_RATIO)>0.4 || min(sp,SPARSITY_RATIO)>0.1)
    if ratio>10
        gamma = 2;
    elseif ratio>1 && ratio<=10
        gamma = 1;
    elseif ratio>0.5 && ratio<=1
        gamma = 0.9;
    elseif ratio>0.1 && ratio<=0.5
        gamma = 0.8;
    elseif ratio>0.05 && ratio<=0.1
        gamma = 0.7;
    elseif ratio>0.01 && ratio<=0.05
        gamma = 0.6;
    elseif ratio>0.005 && ratio<=0.01
        gamma = 0.5;
    elseif ratio>0.001 && ratio<=0.005
        gamma = 0.4;
    else
        gamma = 0.3;
    end
    
    if error_ratio>6 && error_ratio<=10
        sigma = 0.5;  
    elseif error_ratio>3 && error_ratio<=6
        sigma = 0.6; 
        gamma=gamma*0.8;           %add
    elseif error_ratio>1 && error_ratio<=3 && ratio > 0.1
        gamma = gamma*0.4;
        final_check = 1;
        if ratio<0.1
            sigma=0.8;
        else
            sigma=0.7;
        end
    end
else
    if ratio>10
        gamma = 3;
    elseif ratio>1 && ratio<=10
        gamma = 1;
    elseif ratio>0.5 && ratio<=1
        gamma = 0.9;
    elseif ratio>0.1 && ratio<=0.5
        gamma = 0.8;
    elseif ratio>0.05 && ratio<=0.1
        gamma = 0.7;
    elseif ratio>0.01 && ratio<=0.05
        gamma = 0.6;
    elseif ratio>0.005 && ratio<=0.01
        gamma = 0.5;
    elseif ratio>0.001 && ratio<=0.005
        gamma = 0.4;
    else
        gamma = 0.3;
    end    
    
     if error_ratio>6 && error_ratio<=10
        sigma = 0.82;
        gamma=gamma*0.8;
    elseif error_ratio>4 && error_ratio<=6
        sigma = 0.84; 
        gamma=gamma*0.6;           
    elseif error_ratio>3 && error_ratio<=4
        sigma = 0.85; 
        gamma=gamma*0.5;     
        final_check=1;
    elseif error_ratio>1 && error_ratio<=3 
        final_check=1;
        if ratio<0.1
            if double_check
                sigma=0.9;
                gamma=gamma*0.4;
                double_check=0;
            else
                sigma=1.0;
                gamma=gamma*0.1;
                double_check=1;
            end
        else
            sigma=0.88;
            gamma=gamma*0.4;
        end
     else
         sigma=sigma;
     end
end

sigma = sigma*0.8;
mu=mu*sigma;                     
     



      barrier.sigma=sigma;
      barrier.gamma=gamma;
      barrier.mu=mu;
      barrier.final_check=final_check;
      barrier.double_check=double_check;
      barrier.ratio=ratio;

    
end