function [d,val,lam,k]=trustq(gk,Bk,dta)
% 功能: 求解信赖域子问题: min qk(d)=gk'*d+0.5*d'*Bk*d, s.t.||d||<=delta 
%输入: gk是xk处的梯度, Bk是第k次近似Hesse阵, dta是当前信赖域半径
%输出: d, val分别是子问题的最优点和最优值, lam是乘子值, k是迭代次数. 
n=length(gk); 
gamma=0.05;
epsilon=1.0e-6;  
rho=0.6;  
sigma=0.2;
mu0=0.05;  
lam0=0.05;
d0=ones(n,1);  
u0=[mu0,zeros(1,n+1)]';
z0=[mu0,lam0,d0']';
k=0; %k为迭代次数
z=z0; mu=mu0; lam=lam0; d=d0; 
while (k<=150)
    dh=dah(mu,lam,d,gk,Bk,dta);
    if(norm(dh)<epsilon)
        break; 
    end
     A=JacobiH(mu,lam,d,Bk,dta);
     b=beta(mu,lam,d,gk,Bk,dta,gamma)*u0-dh;
     B=inv(A);   
     dz=B*b;
     dmu=dz(1); dlam=dz(2); dd=dz(3:n+2);
     m=0;  mk=0;
     while (m<20)
         dhnew=dah(mu+rho^m*dmu,lam+rho^m*dlam,d+rho^m*dd,gk,Bk,dta);
         if(norm(dhnew)<=(1-sigma*(1-gamma*mu0)*rho^m)*dh)
             mk=m;
             break;
         end
         m=m+1; 
     end
     alpha=rho^mk;
     mu=mu+alpha*dmu;
     lam=lam+alpha*dlam;
     d=d+alpha*dd;
     k=k+1;
end
val=gk'*d+0.5*d'*Bk*d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=phi(mu,a,b)
p=a+b-sqrt((a-b)^2+4*mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dh=dah(mu,lam,d,gk,Bk,dta)
n=length(d);
dh(1)=mu;  dh(2)=phi(mu,lam, dta^2-norm(d)^2);
mh=(Bk+lam*eye(n))*d+gk;
for(i=1:n)
    dh(2+i)=mh(i);
end
dh=dh(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bet=beta(mu,lam,d,gk,Bk,dta,gamma)
dh=dah(mu,lam,d,gk,Bk,dta);
bet=gamma*norm(dh)*min(1,norm(dh));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=JacobiH(mu,lam,d,Bk,dta)
n=length(d);
A=zeros(n+2,n+2);
pmu=-4*mu/sqrt((lam+norm(d)^2-dta^2)^2+4*mu^2);
thetak=(lam+norm(d)^2-dta^2)/sqrt((lam+norm(d)^2-dta^2)^2+4*mu^2);
A=[1,0,zeros(1,n);pmu,1-thetak,-2*(1+thetak)*d'; zeros(n,1),d,Bk+lam*eye(n)];

