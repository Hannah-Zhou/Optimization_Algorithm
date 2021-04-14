function [x,val,k]=sr1(fun,gfun, x0)
%����: �öԳ���1�㷨�����Լ������: min f(x)
%����: x0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�������ݶ� 
%���: x, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������. 
maxk=500; %��������������
rho=0.55;sigma=0.4; epsilon=1e-5;
k=0; n=length(x0); Hk=eye(n);
while(k<maxk)
    gk=feval(gfun,x0); %�����ݶ�
    dk=-Hk*gk; %������������
    if(norm(gk)<epsilon), break; end %������ֹ׼�� 
    m=0; mk=0;
    while(m<20) % ��Armijo�����󲽳�
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
                mk=m; break;
        end
        m=m+1;
    end
    x=x0+rho^mk*dk;
    sk=x-x0; 
    yk=feval(gfun,x)-gk; 
    Hk=Hk+(sk-Hk*yk)*(sk-Hk*yk)'/((sk-Hk*yk)'*yk); %��1У�� 
    k=k+1; x0=x;
end
val=feval(fun,x0);

    