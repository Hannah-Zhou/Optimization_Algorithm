function [x,mu,lam,val,k]=sqpm(x0,mu0,lam0)
%功能: 用基于拉格朗日函数Hesse阵的SQP方法求解约束优化问题: 
% min f(x) s.t. h ?i(x)=0, i=1,..., l.
%输入: x0是初始点, mu0是乘子向量的初始值
%输出: x, mu分别是近似最优点及相应的乘子,
%val是最优值, mh是约束函数的模, k是迭代次数.
maxk=100; %最大迭代次数
n=length(x0); l=length(mu0); m=length(lam0);
rho=0.5; eta=0.1; B0=eye(n);
x=x0; mu=mu0; lam=lam0;
Bk=B0; sigma=0.8;
epsilon1=1e-6; epsilon2=1e-5;
[hk,gk]=cons(x); dfk=df1(x);
[Ae,Ai]=dcons(x); Ak=[Ae; Ai];
k=0;
while(k<maxk)
    [dk,mu,lam]=qpsubp(dfk,Bk,Ae,hk,Ai,gk); %求解子问题 
    mp1=norm(hk,1)+norm(max(-gk,0),1);
    if(norm(dk,1)<epsilon1)&(mp1<epsilon2)
        break;
    end %检验终止准则
    deta=0.05; %罚参数更新 
    tau=max(norm(mu,inf),norm(lam,inf)); 
    if(sigma*(tau+deta)?1)
    sigma=sigma;
    else
        sigma=1.0/(tau+2*deta);
    end
    im=0; %Armijo搜索 
    while(im<=20)
        if(phi1(x+rho^im*dk,sigma)-phi1(x,sigma)<eta*rho^im*dphi1(x,sigma,dk))
            mk=im;
            break; 
        end
        im=im+1;
        if(im==20),mk=10; end
    end
    alpha=rho^mk; 
    x1=x+alpha*dk;
    [hk,gk]=cons(x1);
    dfk=df1(x1);
    [Ae,Ai]=dcons(x1); 
    Ak=[Ae; Ai]; 
    lamu=pinv(Ak)'*dfk; %计算最小二乘乘子
    if(l>0&m>0)
        mu=lamu(1:l); lam=lamu(l+1:l+m);
    end
    if(l==0), mu=[]; lam=lamu; end 
    if(m==0), mu=lamu; lam=[]; end 
    sk=alpha*dk; %更新矩阵Bk 
    yk=dlax(x1,mu,lam)-dlax(x,mu,lam); 
    if(sk'*yk>0.2*sk'*Bk*sk)
        theta=1;
    else
        theta=0.8*sk'*Bk*sk/(sk'*Bk*sk-sk'*yk);
    end
    zk=theta*yk+(1-theta)*Bk*sk;
    Bk=Bk+zk*zk'/(sk'*zk)-(Bk*sk)*(Bk*sk)'/(sk'*Bk*sk);
    x=x1;  k=k+1;
end
val=f1(x);
%p=phi1(x,sigma)
%dd=norm(dk)
%%%%%%%%%%%%%%%%%%%%%% l1精确价值函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=phi1(x,sigma)
f=f1(x); [h,g]=cons(x); gn=max(-g,0); 
l0=length(h); m0=length(g);
if(l0==0), p=f+1.0/sigma*norm(gn,1); end 
if(m0==0), p=f+1.0/sigma*norm(h,1); end 
if(l0>0&m0>0)
    p=f+1.0/sigma*(norm(h,1)+norm(gn,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 价值函数的方向导数%%%%%%%%%%%%%%%%%%%%%%%%%%
function dp=dphi1(x,sigma,d)
df=df1(x); [h,g]=cons(x); gn=max(-g,0); 
l0=length(h); m0=length(g);
if(l0==0), dp=df'*d-1.0/sigma*norm(gn,1); end 
if(m0==0), dp=df'*d-1.0/sigma*norm(h,1); end 
if(l0>0&m0>0)
    dp=df'*d-1.0/sigma*(norm(h,1)+norm(gn,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%% 拉格朗日函数 L(x,mu) %%%%%%%%%%%%%%%%%%%%%%%%%% 
function l=la(x,mu,lam)
f=f1(x); %调用目标函数文件
[h,g]=cons(x); %调用约束函数文件 l0=lemgth(h); m0=length(g);
if(l0==0), l=f-lam*g;  end
if(m0==0),  l=f-mu'*h;  end
if(l0>0&m0>0)
    l=f-mu'*h-lam'*g;
end
%%%%%%%%% 拉格朗日函数的梯度 %%%%%%%%%%%%% 
function dl=dlax(x,mu,lam)
df=df1(x); %调用目标函数梯度文件 
[Ae,Ai]=dcons(x); %调用约束函数Jacobi矩阵文件
[m1,m2]=size(Ai); [l1,l2]=size(Ae);
if(l1==0), dl=df-Ai'*lam; end
if(m1==0), dl=df-Ae'*mu; end
if(l1>0&m1>0), dl=df-Ae>*mu-Ai>*lam; end

