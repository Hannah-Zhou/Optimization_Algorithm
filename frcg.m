function [x,val,k]=frcg(fun,gfun,x0)
% 功能: 用FR共轭梯度法求解无约束问题: min f(x) 
%输入: x0是初始点, fun, gfun分别是目标函数和梯度 
%输出: x, val分别是近似最优点和最优值, k是迭代次数. 
maxk=5000; %最大迭代次数
rho=0.6;sigma=0.4;
k=0; epsilon=1e-4;
n=length(x0);
while(k<maxk)
    g=feval(gfun,x0); %计算梯度
    itern=k-(n+1)*floor(k/(n+1));
    itern=itern+1;
    %计算搜索方向
    if(itern==1)
        d=-g; 
    else
        beta=(g'*g)/(g0'*g0);
        d=-g+beta*d0;  
        gd=g'*d;
        if(gd>=0.0)
            d=-g; 
        end
    end
    if(norm(g)<epsilon)
        break; 
    end %检验终止条件
    m=0; mk=0;
    while(m<20) %Armijo搜索
        if(feval(fun,x0+rho^m*d)<feval(fun,x0)+sigma*rho^m*g'*d)
                 mk=m; break;
        end
        m=m+1; 
    end
    x0=x0+rho^mk*d;
    val=feval(fun,x0);
    g0=g;  
    d0=d;
    k=k+1;
end
x=x0;
val=feval(fun,x);
    