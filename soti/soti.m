function [S,u,v]=soti(Lx,Ly)

xbc=1;
ybc=1;
m=1;
lam=1;
del=0.75;
t=1;

z0=zeros(2);s0=eye(2);s1=[0,1;1,0];s2=[0,-1i;1i,0];s3=[1,0;0,-1]; %#ok<*NASGU>
tzs0=[s0,z0;z0,-s0];
txsx=[z0,s1;s1,z0];
txsy=[z0,s2;s2,z0];
tys0=[z0,-1i*s0;1i*s0,z0];

%Lx=32;
%Ly=24;
N=Lx*Ly*4;
h=zeros(N);

for x=1:Lx
    for y=1:Ly
        i=(x-1)*4+(y-1)*4*Lx;
        
        h(i+1:i+4,i+1:i+4)=m*tzs0;
        
        if x<Lx
            j=x*4+(y-1)*4*Lx;
            boundary=1;
        else
            j=(y-1)*4*Lx;
            boundary=xbc;
        end
        h(i+1:i+4,j+1:j+4)=(0.5*t*tzs0-0.5i*lam*txsx+0.5*del*tys0)*boundary;
        h(j+1:j+4,i+1:i+4)=h(i+1:i+4,j+1:j+4)';
        
        if y<Ly
            j=(x-1)*4+y*4*Lx;
            boundary=1;
        else
            j=(x-1)*4;
            boundary=ybc;
        end
        h(i+1:i+4,j+1:j+4)=(0.5*t*tzs0-0.5i*lam*txsy-0.5*del*tys0)*boundary;
        h(j+1:j+4,i+1:i+4)=h(i+1:i+4,j+1:j+4)';
        
    end
end
[u,v]=eig(h);
%plot(diag(v),'.');


f=diag([ones(1,N/2),zeros(1,N/2)]);
C=u*f*u';
C=(C+C')/2;
%C=zeros(N);
%for i=1:N
%    for j=1:N
%        C(i,j)=u(j,1:N/2)*u(i,1:N/2)';
%    end
%end

sub1=1:N/2;
sub2=[];
sub3=[];
for y=1:Ly
    sub2=[sub2,(1:Lx*2)+(y-1)*Lx*4]; %#ok<*AGROW>
    if y<=Ly/2
        sub3=[sub3,(1:Lx*2)+(y-1)*Lx*4];
    else
        sub3=[sub3,(Lx*2+1:Lx*4)+(y-1)*Lx*4];
    end
end

f1=eig(C(sub1,sub1));
f2=eig(C(sub2,sub2));
f3=eig(C(sub3,sub3));

S = -sum(f1.*log(f1+1e-10)+(1-f1).*log(1-f1+1e-10)) ...
    -sum(f2.*log(f2+1e-10)+(1-f2).*log(1-f2+1e-10))...
    +sum(f3.*log(f3+1e-10)+(1-f3).*log(1-f3+1e-10));
end
        