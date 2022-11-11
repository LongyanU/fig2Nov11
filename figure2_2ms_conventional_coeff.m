clear;
clc
close all

global v ratio M h tau

v=5000;

h=15;
tau=0.002;
M=3;
r=v*tau/h;
ratio=0.8;

x0=0.01*ones(1,M+2);%ϵ���ĳ�ֵ��0
% x0=[ 0.55651, 0.159946, -0.00839528, 0.00251221, 0.201776]

options = optimset('Algorithm','trust-region-reflective','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',2000,'MaxIter',2000);
% options = optimset('Algorithm','interior-point','TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',2000,'MaxIter',2000);

aa=-2*ones(1,M+1);
bb=2*ones(1,M+1);

tic
% [x,resnorm] = lsqnonlin(@myfun7,x0,aa,bb,options)   % Invoke optimizer
[x,resnorm] = lsqnonlin(@myfun7,x0)   % Invoke optimizer
% [x,resnorm] = fmincon(@myfun7,x0,[],[],[],[],aa,bb,[],options)   % Invoke optimizer
toc

% x=real(x)+imag(x);
k=linspace(1/10000,pi/h,50); %

F=zeros(5,50);

for i=1:5  % F����˼��5���Ƕ�,30�������㣬����5���Ƕ�,30�������� ����ʱ��-�ռ���Ƶɢ��ϵ��ʵ�ʾ�������Щ�ض��ĽǶȺͲ����㣬��ʽ�������С��
    xita=(i-1)*pi/16;
    
    for m=1:M
        F(i,:)=2*x(m)*( sin((m-0.5)*k*h *cos(xita) ))+ F(i,:);
    end
    F(i,:)=F(i,:)+4*x(M+1)*cos(k*sin(xita)*h).*( sin(k*h*cos(xita)/2) );
    F(i,:)=F(i,:).*(1+2*x(M+2)*(cos(k*h*sin(xita))-1));
    F(i,:)= F(i,:).^2;
    
    temp=0;
    for m=1:M
        temp=2*x(m)*( sin( (m-0.5)*k*h *sin(xita) ) )+ temp;
    end
    temp=temp+4*x(M+1)*cos(k*cos(xita)*h).*( sin(k*h*sin(xita)/2) );
    temp=temp.*(1+2*x(M+2)*(cos(k*h*cos(xita))-1));
    
    F(i,:)=F(i,:)+temp.^2;
    F(i,:)=F(i,:)*r^2;
    

    F(i,:)=-1/2*F(i,:)./  ( (1+2*x(M+2)*(cos(k*h*cos(xita))-1)).^2.*(1+2*x(M+2)*(cos(k*h*sin(xita))-1)).^2 )  +1;
    F(i,:)=(acos(F(i,:))./ (tau*k*v));
    a1=(h/v*(1./F(i,:)-1));
    if (i==1)
        figure;plot(v*k*h/(2*pi*h),a1,'m','linewidth',2)
        hold on
    elseif i==2
        plot(v*k*h/(2*pi*h),a1,'r--','linewidth',2)
    elseif i==3
        plot(v*k*h/(2*pi*h),a1,'c:','linewidth',2)
    elseif i==4
        plot(v*k*h/(2*pi*h),a1,'k-.','linewidth',2)
    else
        plot(v*k*h/(2*pi*h),a1,'b','linewidth',2)
    end
end
grid on

% axis([0 v/(2*h) -3*10^-5 7*10^-5])
xlabel('f(hz)')
legend('\theta=0', '\theta=��/16','\theta=2��/16','\theta=3��/16','\theta=4��/16')

ylabel('\delta (\theta)')