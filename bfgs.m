function [x,val,k]=bfgs(fun,gfun,x0,varargin)
%����: ��BFGS�㷨�����Լ������: min f(x)
%����: x0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�������ݶ�;
% varargin������Ŀɱ��������, �򵥵���bfgsʱ���Ժ�����, 
% ������������ѭ�����øó���ʱ��������Ҫ������
%���: x, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������. 
maxk=500; %��������������
rho=0.55;sigma=0.4; epsilon=1e-5;
k=0; n=length(x0);
 Bk=eye(n);  %Bk=feval(��Hess��,x0);
 while(k<maxk)
     gk=feval(gfun,x0,varargin{:}); %�����ݶ� 
     if(norm(gk)<epsilon), break; end %������ֹ׼�� 
     dk=-Bk\gk; %�ⷽ����, ������������
     m=0; mk=0;
     while(m<20) % ��Armijo�����󲽳� 
         newf=feval(fun,x0+rho^m*dk,varargin{:}); 
         oldf=feval(fun,x0,varargin{:}); 
         if(newf<oldf+sigma*rho^m*gk'*dk)
             mk=m; break;
         end
         m=m+1;
     end
     %BFGSУ��
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