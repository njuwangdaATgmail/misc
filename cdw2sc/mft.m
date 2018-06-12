function[delta,mu]=mft(T,L,V1,V2,Vnn,doping)
% calling from outside
% for i=1:length(doping)
% for j=1:nsample
% [delta,mu]=mft(0.05, 16, 2, 0, 0.03, doping(i));
% sc_max(i,j)=max(max(abs(delta)));
% sc_mean(i,j)=mean(mean(abs(delta)));
% cdw_max(i,j)=max(max(abs(mu)));
% cdw_mean(i,j)=mean(mean(abs(mu)));
% disp([i,j])
% end
% end

t=1;
mu0=0.1;
converg=1e-4;
mixing=0.5;
iter_max=20000;

N=L*L;
M=floor(N*doping);
Vph=V1*ones(1,N);
if M>0
    Vph(randperm(N,M))=V2;
end

h0=zeros(N);
for x=1:L
    for y=1:L
        i=x+(y-1)*L;
        h0(i,i)=-mu0;
        xx=mod(x,L)+1;
        yy=y;
        ii=xx+(yy-1)*L;
        h0(i,ii)=-t;
        h0(ii,i)=h0(i,ii)';
        xx=x;
        yy=mod(y,L)+1;
        ii=xx+(yy-1)*L;
        h0(i,ii)=-t;
        h0(ii,i)=h0(i,ii)';
    end
end

mu=rand(1,N)-0.5;
delta=rand(1,N)-0.5;
err=converg+1;
iter=0;
while err>converg
    
    iter=iter+1;
    if iter>iter_max
        disp('not converged');
        break
    end
    
    H=[h0+diag(mu),diag(delta);
        diag(delta'),-h0-diag(mu)];
    [U,V]=eig(H);
    th=diag(tanh(diag(V(1:N,1:N))/2/T));
    u=U(1:N,1:N);
    v=U(N+1:2*N,1:N);
    deltanew=diag(u*th*v')'.*Vph;
    ne=(diag(u*th*u'-v*th*v')'); %.*Vph/2; %+ones(1,N)
    ne=ne-mean(ne);
    munew=ne.*Vph/2;
    for x=1:L
        for y=1:L
            i=x+(y-1)*L;
            xx=mod(x,L)+1;
            yy=y;
            ii=xx+(yy-1)*L;
            munew(i)=munew(i)-Vnn/2*ne(ii);
            xx=mod(x-2,L)+1;
            yy=y;
            ii=xx+(yy-1)*L;
            munew(i)=munew(i)-Vnn/2*ne(ii);
            xx=x;
            yy=mod(y,L)+1;
            ii=xx+(yy-1)*L;
            munew(i)=munew(i)-Vnn/2*ne(ii);
            xx=x;
            yy=mod(y-2,L)+1;
            ii=xx+(yy-1)*L;
            munew(i)=munew(i)-Vnn/2*ne(ii);
        end
    end
    %munew(1:4)
    %pause
    
    err=max(abs([deltanew-delta,munew-mu]));
    if mod(iter,10000)==0
        disp([iter,err]);
    end
    
    delta=delta+(deltanew-delta)*mixing;
    mu=mu+(munew-mu)*mixing;
    
end 
delta=reshape(delta,L,L);
mu=reshape(mu,L,L)-mean(mean(mu));
end
