function [s,phis,k,G,E]=golds(phi,a,b,delta,epsilon) 
%����: phi��Ŀ�꺯��, a, b ����������������˵�
% delta, epsilon�ֱ����Ա����ͺ���ֵ��������� %���: s, phis�ֱ��ǽ��Ƽ�С��ͼ�Сֵ, G��nx4����,
% ���k�зֱ���a,p,q,b�ĵ�k�ε���ֵ[ak,pk,qk,bk],
% E=[ds,dphi], �ֱ���s��phis�������. t=(sqrt(5)-1)/2; h=b-a; phia=feval(phi,a); phib=feval(phi,b);
p=a+(1-t)*h;  
q=a+t*h;
phip=feval(phi,p); 
phiq=feval(phi,q);
k=1;  
G(k,:)=[a, p, q, b];
while(abs(phib-phia)>epsilon)|(h>delta)
    if(phip<phiq)
        b=q;  
        phib=phiq; 
        q=p; 
        phiq=phip;
        h=b-a; 
        p=a+(1-t)*h; 
        phip=feval(phi,p);
    else
        a=p; 
        phia=phip; 
        p=q; 
        phip=phiq;
        h=b-a;  
        q=a+t*h;  
        phiq=feval(phi,q);
end
    k=k+1;  
    G(k,:)=[a, p, q, b];
end
ds=abs(b-a); dphi=abs(phib-phia);
if(phip<=phiq)
    s=p;  phis=phip;
else
    s=q;  phis=phiq;
end
E=[ds,dphi];