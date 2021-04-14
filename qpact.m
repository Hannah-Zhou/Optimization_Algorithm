function [x,lamk,exitflag,output]=qpact(H,c,Ae,be,Ai,bi,x0)
%功能: 用有效集方法解一般约束二次规划问题:
% min f(x)=0.5*x'*H*x+c'*x,
% s.t. a'_i*x-b_i=0,(i=1,...,l),
%      a'_i*x-b_i>=0,(i=l+1,...,m)
%输入: x0是初始点, H, c分别是目标函数二次型矩阵和向量;
% Ae=(a_1,...,a_l)’, be=(b_1,...,b_l)';
% Ai=(a_{l+1},...,a_m), bi=(b_{l+1},...,n_m)'.
%输出: x是最优解,lambda是对应的乘子向量;output是结构变量,输出极小值f(x),迭代次数k等信息,exitflag是算法终止类型
%%%%%%%%%%%%%%%%% 主程序开始 %%%%%%%%%%%%%%%%%
% 初始化
epsilon=1.0e-9; err=1.0e-6;
k=0; x=x0; n=length(x); kmax=1.0e3; 
ne=length(be); ni=length(bi); lamk=zeros(ne+ni,1); 
index=ones(ni,1);
for (i=1:ni)
    if(Ai(i,:)*x>bi(i)+epsilon), index(i)=0; end
end
%算法主程序 
while (k<=kmax)
    %求解子问题
    Aee=[];
    if(ne>0), Aee=Ae; end 
    for(j=1:ni)
        if(index(j)>0), Aee=[Aee; Ai(j,:)]; end
    end
    gk=H*x+c;
    [m1,n1] = size(Aee);
    [dk,lamk]=qsubp(H,gk,Aee,zeros(m1,1));
    if(norm(dk)<=err)
        y=0.0;
        if(length(lamk)>ne)
            [y,jk]=min(lamk(ne+1:length(lamk)));
        end
        if(y>=0)
            exitflag=0;
        else
            exitflag=1;
            for(i=1:ni)
                if(index(i)&(ne+sum(index(1:i)))==jk)
                    index(i)=0; break;
                end
            end
        end
        k=k+1; 
    else
        exitflag=1; 
        %求步长
        alpha=1.0; tm=1.0; 
        for(i=1:ni)
                 if((index(i)==0)&(Ai(i,:)*dk<0))
                     tm1=(bi(i)-Ai(i,:)*x)/(Ai(i,:)*dk);
                     if(tm1<tm)
                         tm=tm1; ti=i;
                     end
                 end
        end
        alpha=min(alpha,tm); 
        x=x+alpha*dk;
        %修正有效集
        if(tm<1), index(ti)=1; end
    end
         if(exitflag==0), break; end
         k=k+1;
end
output.fval=0.5*x'*H*x+c'*x; 
output.iter=k;
%%%%%%%%%%%%%%%%%%%%%%%%%% 求解子问题 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [x,lambda]=qsubp(H,c,Ae,be) 
ginvH=pinv(H);
[m,n]=size(Ae);
if(m>0)
    rb=Ae*ginvH*c + be;
    lambda=pinv(Ae*ginvH*Ae')*rb;
    x=ginvH*(Ae'*lambda-c);
else
    x=-ginvH*c;
    lambda=0;
end

