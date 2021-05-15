function vars=reinitialize_vars(sigma, vars, index)

m=vars.m;
l=vars.l;
u=vars.u;
v=vars.v;

if index==0
    for i=m+1:l
        if u(i)>v(i)
            v(i)=sigma*v(i);
        else
            u(i)=sigma*u(i);
        end
    end
elseif index==1
    for i=m+1:l
        u(i)=sqrt(sigma)*u(i);
        v(i)=sqrt(sigma)*v(i);
    end
else
    for i=m+1:l
        u(i)=sqrt(1/sigma)*u(i);
        v(i)=sqrt(1/sigma)*v(i);
    end    
end

vars.u=u;
vars.v=v;


end
    
    
    















