function [xk,val,k]=trustm(x0)
%����: ţ���������򷽷������Լ���Ż����� min f(x) 
%����: x0�ǳ�ʼ������
%���: xk�ǽ��Ƽ�С��, val�ǽ��Ƽ�Сֵ, k�ǵ������� 
n=length(x0); x=x0; dta=1;
eta1=0.1; eta2=0.75; dtabar=2.0;
tau1=0.5; tau2=2.0; epsilon=1e-6;
k=0; Bk=Hess(x);
while(k<50)
         gk=gfun(x);
         if(norm(gk)<epsilon)
             break; 
         end
         % �����ӳ���trustq 
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

