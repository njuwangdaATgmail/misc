function[G0,Sig0,G0exact,Sig0exact,closed]=gaussianphi(beta,Nk)
g=1;

dphi=0.001;
phi=0:dphi:2*max(1/beta,g);  % |phi|^2>>T*g, where g=1
iw0=1i*pi/beta;

G11=0.5*(1./(iw0-phi)+1./(iw0+phi));
G12=0.5*(1./(iw0-phi)-1./(iw0+phi));

for iphi=1:length(phi)
    if beta*phi(iphi)<10
        S0(iphi)=2*log(2+2*cosh(beta*phi(iphi)))-beta/g*phi(iphi)^2;
    else
        S0(iphi)=2*beta*phi(iphi)-beta/g*phi(iphi)^2;
    end
end
%S0=2*log(2+2*cosh(beta*phi))-beta/g*phi.^2;

if Nk==0
    closed=1;
    Gphi=G11;
elseif Nk==1
    disp('incorrect now...');return
    closed=1-0.5*g*beta+0.5*g*beta*tanh(beta*phi/2).^2;
    Gphi=G11.*closed+G11*0.5*g.*G11+G12*0.5*g.*G12+G11.*G12*g.*tanh(beta*phi/2); 
end

nv=(1:0*g*beta/2/pi)*2*pi/beta;
dS=zeros(1,length(phi));
for iphi=1:length(phi)
    D=1-4*g*phi(iphi)./(nv.^2+4*phi(iphi)^2)*tanh(beta*phi(iphi)/2);
    sgn=sign(D);
    dSnv=-log(abs(D));
    dS(iphi)=sum(dSnv);
    Gphi(iphi)=Gphi(iphi)*exp(sum(log(sgn)));
end

Wphi=phi.*exp(S0+dS-max(S0+dS));
G0=sum(Gphi.*Wphi)/sum(closed.*Wphi);
Sig0=iw0-1/G0;

%plot(phi,Wphi);

Z=14+2*cosh(beta);
G0exact=((1/(iw0+1)+1/(iw0-1))*(1+cosh(beta))+12/iw0)/Z;
Sig0exact=iw0-1/G0exact;

end
