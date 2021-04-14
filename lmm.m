function [x,val,k]=lmm(Fk,JFk,x0)
%����: ��L-M�����������Է�����: F(x)=0
%����: x0�ǳ�ʼ��, Fk, JFk �ֱ�����F(xk)��F��(xk)�ĺ���. 
%���: x, val�ֱ��ǽ��ƽ⼰����F(xk)������ֵ, k�ǵ�������. 
maxk=100; %��������������
rho=0.55;sigma=0.4; muk=norm(feval(Fk,x0));
k=0; epsilon=1e-6; n=length(x0);
while(k<maxk)
    fk=feval(Fk,x0); %���㺯��ֵ
    jfk=feval(JFk,x0); %����Jacobi��
    gk=jfk'*fk;
    dk=-(jfk'*jfk+muk*eye(n))\gk; %�ⷽ����Gk*dk=-gk, ������������ 
    if(norm(gk)<epsilon), break; end %������ֹ׼��
    m=0; mk=0;
    while(m<20) % ��Armijo�����󲽳�
        newf=0.5*norm(feval(Fk,x0+rho^m*dk))^2;
        oldf=0.5*norm(feval(Fk,x0))^2;
        if(newf<oldf+sigma*rho^m*gk;*dk)
            mk=m; 
            break;
        end
        m=m+1; 
    end
    x0=x0+rho^mk*dk;
    muk=norm(feval(Fk,x0));
    k=k+1;
end
x=x0;
val=0.5*muk^2;