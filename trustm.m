function [xk,val,k]=trustm(x0)
%功能: 牛顿型信赖域方法求解无约束优化问题 min f(x) 
%输入: x0是初始迭代点
%输出: xk是近似极小点, val是近似极小值, k是迭代次数 
n=length(x0); x=x0; dta=1;
eta1=0.1; eta2=0.75; dtabar=2.0;
tau1=0.5; tau2=2.0; epsilon=1e-6;
k=0; Bk=Hess(x);
while(k<50)
         gk=gfun(x);
         if(norm(gk)<epsilon)
             break; 
         end
         % 调用子程序trustq 
         [d,val,lam,ik]=trustq(gk,Bk,dta); 
         deltaq=-qk(x,d); 
         deltaf=fun(x)-fun(x+d); 
         rk=deltaf/deltaq;
         if(rk<=eta1)
             dta=tau1*dta;
         else if (rk>=eta2&norm(d)==dta)
                 dta=min(tau2*dta,dtabar);
             else
                 dta=dta; 
             end
         end
         if(rk>eta1)
             x=x+d;
       Bk=Hess(x);
         end
         k=k+1; 
end
xk=x;
val=fun(xk);

