function [x,val,k]=broyden(fun,gfun,x0)
%功能: 用Broyden族算法求解无约束问题: min f(x) 
%输入: x0是初始点, fun, gfun分别是目标函数及其梯度 
%输出: x,val分别是近似最优点和最优值, k是迭代次数. 
maxk=1e5; %给出最大迭代次数
rho=0.55;sigma=0.4; epsilon=1e-5;
phi=0.5; k=0; n=length(x0); Hk=inv(feval('Hess',x0)); %Hk=eye(n); 
while(k<maxk)'
    gk=feval(gfun,x0); %计算梯度 
    if(norm(gk)<epsilon), break; end %检验终止准则 
    dk=-Hk*gk; %解方程组,计算搜索方向
    m=0; mk=0;
    while(m<20) % 用Armijo搜索求步长
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
            mk=m; break;
        end
        m=m+1; 
    end
    %Broyden族校正
    x=x0+rho^mk*dk;
    sk=x-x0; 
    yk=feval(gfun,x)-gk; 
    Hy=Hk*yk; sy=sk'*yk;
    yHy=yk'*Hk*yk; 
    if(sy<0.2*yHy)
        theta=0.8*yHy/(yHy-sy);
        sk=theta*sk+(1-theta)*Hy;
        sy=0.2*yHy;
    end
    vk=sqrt(yHy)*(sk/sy - Hy/yHy);
    Hk=Hk-(Hy*Hy')/yHy+(sk*sk')/sy+phi*vk*vk';
    k=k+1;  
    x0=x;
end
val=feval(fun,x0);

