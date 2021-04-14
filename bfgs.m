function [x,val,k]=bfgs(fun,gfun,x0,varargin)
%功能: 用BFGS算法求解无约束问题: min f(x)
%输入: x0是初始点, fun, gfun分别是目标函数及其梯度;
% varargin是输入的可变参数变量, 简单调用bfgs时可以忽略它, 
% 但若其它程序循环调用该程序时将发挥重要的作用
%输出: x, val分别是近似最优点和最优值, k是迭代次数. 
maxk=500; %给出最大迭代次数
rho=0.55;sigma=0.4; epsilon=1e-5;
k=0; n=length(x0);
 Bk=eye(n);  %Bk=feval(’Hess’,x0);
 while(k<maxk)
     gk=feval(gfun,x0,varargin{:}); %计算梯度 
     if(norm(gk)<epsilon), break; end %检验终止准则 
     dk=-Bk\gk; %解方程组, 计算搜索方向
     m=0; mk=0;
     while(m<20) % 用Armijo搜索求步长 
         newf=feval(fun,x0+rho^m*dk,varargin{:}); 
         oldf=feval(fun,x0,varargin{:}); 
         if(newf<oldf+sigma*rho^m*gk'*dk)
             mk=m; break;
         end
         m=m+1;
     end
     %BFGS校正
     x=x0+rho^mk*dk;
     sk=x-x0; 
     yk=feval(gfun,x,varargin{:})-gk; 
     if(yk'*sk>0)
         Bk=Bk-(Bk*sk*sk'*Bk)/(sk'*Bk*sk)+(yk*yk')/(yk'*sk);
     end
     k=k+1;   
     x0=x;
 end
 val=feval(fun,x0,varargin{:});