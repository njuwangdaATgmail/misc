function [E,y]=soti_ky(Lx,Ly,xbc)

m=1;
lam=1;
del=0.5;
t=1;

z0=zeros(2);s0=eye(2);s1=[0,1;1,0];s2=[0,-1i;1i,0];s3=[1,0;0,-1];
tzs0=[s0,z0;z0,-s0];
txsx=[z0,s1;s1,z0];
txsy=[z0,s2;s2,z0];
tys0=[z0,-1i*s0;1i*s0,z0];

N=Lx*4;
E=zeros(N,Ly);
f1=zeros(N/2,Ly);
for y=1:Ly
    ky=2*pi*y/Ly;
    cosky=cos(ky)
    sinky=sin(ky)

    h=zeros(N);

    for x=1:Lx
        
        i=(x-1)*4;
        
        h(i+1:i+4,i+1:i+4)=m*tzs0+t*cosky*tzs0+lam*sinky*txsy-del*cosky*tys0;
        
        if x<Lx
            j=x*4;
            boundary=1;
        else
            j=0;
            boundary=xbc;
        end
        h(i+1:i+4,j+1:j+4)=(0.5*t*tzs0-0.5i*lam*txsx+0.5*del*tys0)*boundary;
        h(j+1:j+4,i+1:i+4)=h(i+1:i+4,j+1:j+4)';
        
    end
    [u,v]=eig(h);
    E(:,y)=diag(v);

    f=diag([ones(1,N/2),zeros(1,N/2)]);
    C=u*f*u';
    C=(C+C')/2;
    sub=1:N/2;
    f1(:,y)=eig(C(sub,sub));
end

end 
