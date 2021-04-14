function mk=armijo(xk,dk )
    beta=0.5;  
    sigma=0.2;
    m=0; 
    mmax=20;
    while (m<=mmax)
        if(fun(xk+beta^m*dk)<=fun(xk)+sigma*beta^m*gfun(xk)'*dk)
            mk=m;  
            break;
        end
    m=m+1; 
    end
    alpha=beta^mk
    newxk=xk+alpha*dk
    fk=fun(xk)
    newfk=fun(newxk)