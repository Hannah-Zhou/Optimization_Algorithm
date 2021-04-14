function [d,mu,lam,val,k]=qpsubp(dfk,Bk,Ae,hk,Ai,gk)
% 功能: 求解二次规划子问题: min qk(d)=0.5*d'*Bk*d+dfk'*d+,
%                           s.t. hk+Ae*d=0, gk+Ai*d?=0.
%输入: dfk是xk处的梯度, Bk是第k次近似Hesse阵, Ae,hk线性等式约束的有关参数,Ai,gk是线性不等式约束的有关参数
%输出: d,val分别是是最优解和最优值, mu,lam是乘子向量, k是迭代次数. 
n=length(dfk); l=length(hk); m=length(gk);
gamma=0.05; epsilon=1.0e-6; rho=0.5; sigma=0.2;
ep0=0.05; mu0=0.05*zeros(l,1); lam0=0.05*zeros(m,1);
d0=ones(n,1); u0=[ep0;zeros(n+l+m,1)];
z0=[ep0; d0; mu0;lam0,];
k=0; %k为迭代次数
z=z0; ep=ep0; d=d0; mu=mu0; lam=lam0;
while (k<=150)
   dh=dah(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk);
   if(norm(dh)<epsilon)
       break; 
   end
   A=JacobiH(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk);
   b=beta(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk,gamma)*u0-dh;
   dz=A\b;
   if(l>0&m>0)
       de=dz(1); dd=dz(2:n+1); du=dz(n+2:n+l+1); dl=dz(n+l+2:n+l+m+1);
   end
   if(l==0)
      de=dz(1);   dd=dz(2:n+1);  dl=dz(n+2:n+m+1);
   end
   if(m==0)
      de=dz(1);   dd=dz(2:n+1);  du=dz(n+2:n+l+1);
   end
   i=0; %mk=0;
   while (mm<=20)
     if(l>0&m>0)
      dh1=dah(ep+rho^i*de,d+rho^i*dd,mu+rho^i*du,lam+rho^i*dl,dfk,Bk,Ae,hk,Ai,gk);\
     end
     if(l==0)
       dh1=dah(ep+rho^i*de,d+rho^i*dd,mu,lam+rho^i*dl,dfk,Bk,Ae,hk,Ai,gk);
     end
     if(m==0)
   dh1=dah(ep+rho^i*de,d+rho^i*dd,mu+rho^i*du,lam,dfk,Bk,Ae,hk,Ai,gk);
     end
     if(norm(dh1)<=(1-sigma*(1-gamma*ep0)*rho^i)*norm(dh))
         mk=i; break;
     end
     i=i+1;
     if(i==20), mk=10; end
   end
   alpha=rho^mk;
   if(l>0&m>0)
     ep=ep+alpha*de; d=d+alpha*dd;
     mu=mu+alpha*du;lam=lam+alpha*dl;
   end
   if(l==0)
     ep=ep+alpha*de; d=d+alpha*dd;
     lam=lam+alpha*dl;
   end
   if(m==0)
     ep=ep+alpha*de; d=d+alpha*dd;
     mu=mu+alpha*du;
   end
   k=k+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=phi(ep,a,b)
p=a+b-sqrt(a^2+b^2+2*ep^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dh=dah(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk)
n=length(dfk); l=length(hk); m=length(gk);
dh=zeros(n+l+m+1,1);
dh(1)=ep;
if(l>0&m>0)
   dh(2:n+1)=Bk*d-Ae'*mu-Ai'*lam+dfk;
   dh(n+2:n+l+1)=hk+Ae*d;
   for(i=1:m)
       dh(n+l+1+i)=phi(ep,lam(i),gk(i)+Ai(i,:)*d);
   end
end
if(l==0)
   dh(2:n+1)=Bk*d-Ai'*lam+dfk;
   for(i=1:m)
       dh(n+1+i)=phi(ep,lam(i),gk(i)+Ai(i,:)*d);
   end
end
if(m==0)
    dh(2:n+1)=Bk*d-Ae'*mu+dfk;
    dh(n+2:n+l+1)=hk+Ae*d;
end
dh=dh(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bet=beta(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk,gamma)
dh=dah(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk);
bet=gamma*norm(dh)*min(1,norm(dh));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dd1,dd2,v1]=ddv(ep,d,lam,Ai,gk)
m=length(gk);
dd1=zeros(m,m); dd2=zeros(m,m); v1=zeros(m,1);
for(i=1:m)
    fm=sqrt(lam(i)^2+(gk(i)+Ai(i,:)*d)^2+2*ep^2);
    dd1(i,i)=1-lam(i)/fm;
    dd2(i,i)=1-(gk(i)+Ai(i,:)*d)/fm;
    v1(i)=-2*ep/fm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=JacobiH(ep,d,mu,lam,dfk,Bk,Ae,hk,Ai,gk)
n=length(dfk); l=length(hk); m=length(gk);
A=zeros(n+l+m+1,n+l+m+1);
[dd1,dd2,v1]=ddv(ep,d,lam,Ai,gk);
if(l>0&m>0)
    A=[1,            zeros(1,n),  zeros(1,l),  zeros(1,m);
       zeros(n,1),   Bk,          -Ae',        -Ai';
       zeros(l,1),   Ae,         zeros(l,l),   zeros(l,m) ;
       v1,           dd2*Ai,     zeros(m,l),   dd1];
end
if(l==0)
          A=[1,          zeros(1,n),  zeros(1,m);
             zeros(n,1), Bk,          -Ai';
             v1,         dd2*Ai,      dd1];
end
if(m==0)
         A=[1,           zeros(1,n), zeros(1,l);
            zeros(n,1),  Bk,          -Ae';
            zeros(l,1),  Ae,        zeros(l,l)];
end

