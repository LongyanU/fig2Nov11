function F = myfun2(x)

global v ratio M h tau;

r=v*tau/h;

k=linspace(1*pi/(1000),ratio*pi/h,20); %
F=zeros(5,20);

for i=1:5  % F����˼��5���Ƕ�,30�������㣬����5���Ƕ�,30�������� ����ʱ��-�ռ���Ƶɢ��ϵ��ʵ�ʾ�������Щ�ض��ĽǶȺͲ����㣬��ʽ�������С��
    xita=(i-1)*pi/16;
    
    for m=1:M
        F(i,:)=2*x(m)*( cos(m*k*h *cos(xita) )-cos((m-1)*k*h*cos(xita)) )+ F(i,:);
    end
    F(i,:)=F(i,:)+4*x(M+1)*cos(k*sin(xita)*h).*(cos(k*h*cos(xita))-1);
    F(i,:)=F(i,:).*(1+2*x(M+2)*(cos(k*h*sin(xita))-1));
    
    
    temp=0;
    for m=1:M
        temp=2*x(m)*( cos(m*k*h *sin(xita) )-cos((m-1)*k*h*sin(xita)) )+ temp;
    end
    temp=temp+4*x(M+1)*cos(k*cos(xita)*h).*(cos(k*h*sin(xita))-1); 
    temp=temp.*(1+2*x(M+2)*(cos(k*h*cos(xita))-1));
    
    F(i,:)=F(i,:)+temp;
    F(i,:)=F(i,:)*r^2;
    
    F(i,:)=1/2*F(i,:)./  ( (1+2*x(M+2)*(cos(k*h*cos(xita))-1)) .*(1+2*x(M+2)*(cos(k*h*sin(xita))-1)) )  +1;
    F(i,:)=(acos(F(i,:))./ (tau*k*v));
end
F=F-1;