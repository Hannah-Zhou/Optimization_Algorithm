function [x,val,k]=sr1(fun,gfun, x0)
%功能: 用对称秩1算法求解无约束问题: min f(x)
%输入: x0是初始点, fun, gfun分别是目标函数及其梯度 
%输出: x, val分别是近似最优点和最优值, k是迭代次数. 
maxk=500; %给出最大迭代次数
rho=0.55;sigma=0.4; epsilon=1e-5;
k=0; n=length(x0); Hk=eye(n);
while(k<maxk)
    gk=feval(gfun,x0); %计算梯度
    dk=-Hk*gk; %计算搜索方向
    if(norm(gk)<epsilon), break; end %检验终止准则 
    m=0; mk=0;
    while(m<20) % 用Armijo搜索求步长
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
                mk=m; break;
        end
        m=m+1;
    end
    x=x0+rho^mk*dk;
    sk=x-x0; 
    yk=feval(gfun,x)-gk; 
    Hk=Hk+(sk-Hk*yk)*(sk-Hk*yk)'/((sk-Hk*yk)'*yk); %秩1校正 
    k=k+1; x0=x;
end
val=feval(fun,x0);

    