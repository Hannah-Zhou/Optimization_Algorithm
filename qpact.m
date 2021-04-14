function [x,lamk,exitflag,output]=qpact(H,c,Ae,be,Ai,bi,x0)
%����: ����Ч��������һ��Լ�����ι滮����:
% min f(x)=0.5*x'*H*x+c'*x,
% s.t. a'_i*x-b_i=0,(i=1,...,l),
%      a'_i*x-b_i>=0,(i=l+1,...,m)
%����: x0�ǳ�ʼ��, H, c�ֱ���Ŀ�꺯�������;��������;
% Ae=(a_1,...,a_l)��, be=(b_1,...,b_l)';
% Ai=(a_{l+1},...,a_m), bi=(b_{l+1},...,n_m)'.
%���: x�����Ž�,lambda�Ƕ�Ӧ�ĳ�������;output�ǽṹ����,�����Сֵf(x),��������k����Ϣ,exitflag���㷨��ֹ����
%%%%%%%%%%%%%%%%% ������ʼ %%%%%%%%%%%%%%%%%
% ��ʼ��
epsilon=1.0e-9; err=1.0e-6;
k=0; x=x0; n=length(x); kmax=1.0e3; 
ne=length(be); ni=length(bi); lamk=zeros(ne+ni,1); 
index=ones(ni,1);
for (i=1:ni)
    if(Ai(i,:)*x>bi(i)+epsilon), index(i)=0; end
end
%�㷨������ 
while (k<=kmax)
    %���������
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
        %�󲽳�
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
        %������Ч��
        if(tm<1), index(ti)=1; end
    end
         if(exitflag==0), break; end
         k=k+1;
end
output.fval=0.5*x'*H*x+c'*x; 
output.iter=k;
%%%%%%%%%%%%%%%%%%%%%%%%%% ��������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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

