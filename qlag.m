function [x,lam,fval]=qlag(H,A,b,c)
% ����: ���������շ�������ʽԼ�����ι滮:
% min f(x)=0.5*x��Hx+c��x, s.t. Ax=b
%����: H,c�ֱ���Ŀ�꺯���ľ��������, A,b�ֱ���Լ�������еľ��������
%���: (x, lam) �� KT ��, fval ������ֵ. 
IH=inv(H);
AHA=A*IH*A';
IAHA=inv(AHA);
AIH=A*IH;
G=IH-AIH'*IAHA*AIH;
B=IAHA*AIH;
C=-IAHA;
x=B'*b-G*c;
lam=B*c-C*b;
fval=0.5*x'*H*x+c'*x;