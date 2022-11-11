clear
clc
close all

% 历时 199.935852 秒。
nt = 1103;

isnap=10;    % snapshot sampling
nx=900;
nz=900;

v=ones(nz,nx)*5000;

dx=15;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.004; % calculate time step from stability criterion
tau=dt;

f0=35;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^10*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=nz/2;
xs=nx/2;


seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;

r=v*dt/h;
coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
        
        coeff(ii,jj,:)= ...
            [0.692236, 0.0347923, -0.0428603, 0.214669, -0.277533];
% [0.732494, 0.0214622, -0.0453531, 0.223779, -0.309616];

        
        
        
    end
end



b=1-2*coeff(:,:,end);
a = coeff(:,:,end);
cc = a;

for j=1:nx,
    A1{j} = (gallery('tridiag',a(1:nz-1,j),b(:,j),cc(1:nz-1,j)));
end



A2=cell(nz,1);
for k=1:nz,
    A2 {k}= (gallery('tridiag',a(k,1:nx-1),b(k,:),cc(k,1:nx-1) ));      %%解三对角阵，直接matlab解了
end


Vx=zeros(nz,nx);
Vz=zeros(nz,nx);
d2pzz=p;d2pxx=p;
tic
for it=1:nt-2,
    
    d2px11=(circshift(Vx,[0 0])-circshift(Vx,[0 1]));
    d2px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
    d2px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
    
    addPoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]));
    addPoint=addPoint+(circshift(Vx,[1 0])-circshift(Vx,[1 1]));
    d2px=coeff(:,:,1).*d2px11+coeff(:,:,2).*d2px12+coeff(:,:,3).*d2px13;
    d2px=d2px+addPoint.*coeff(:,:,4);
    
    d2pz11=(circshift(Vz,[0 0])-circshift(Vz,[-1 0]));
    d2pz12=(circshift(Vz,[1 0])-circshift(Vz,[-2 0]));
    d2pz13=(circshift(Vz,[2 0])-circshift(Vz,[-3 0]));
    
    addPoint=(circshift(Vz,[0 -1])-circshift(Vz,[-1 -1]));
    addPoint=addPoint+(circshift(Vz,[0 1]) -circshift(Vz,[-1 1]));
    
    d2pz=coeff(:,:,1).*d2pz11+coeff(:,:,2).*d2pz12+coeff(:,:,3).*d2pz13;
    d2pz=d2pz+addPoint.*coeff(:,:,4);
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';
    end
    p=p-dt*v.^2.*(d2pxx+d2pzz)/h;
    p(zs,xs)=p(zs,xs)+(src(it)+src(it+1))/2;
    %     [p,p]=spongeABC(p,p,nx,nz,45,45,0.009);
    % time saved
    d2px1=circshift(p,[0 -1])-circshift(p,[0 0]);
    d2pz1=circshift(p,[1])-p;
    
    Vx=Vx-dt*d2px1/h;
    Vz=Vz-dt*d2pz1/h;
    
%     [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    if rem(it,isnap)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    if it==300
        ptemp= p;
    end
    
    if it==350
        ptemp2= p;
    end
    
    
     if it==500
        ptemp3= p;
     end
    
       if it==700
        ptemp4= p;
    end
    
    %   record1(it)=p(680,680);
    
end


toc
save('figure2_4ms_HEI.mat')