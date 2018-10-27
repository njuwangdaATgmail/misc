function[G0,Sig0,G0exact,Sig0exact]=uniformphi(beta,Nk)
g=1;

dphi=0.01;
phi=0:dphi:100*sqrt(g/beta);  % |phi|^2>>T*g, where g=1
iw0=1i*pi/beta;

G11=0.5*(1./(iw0-phi)+1./(iw0+phi));
G12=0.5*(1./(iw0-phi)-1./(iw0+phi));

if Nk==0
    closed=1;
    Gphi=G11;
elseif Nk==1
    closed=1-0.5*g*beta+0.5*g*beta*tanh(beta*phi/2).^2;
    Gphi=G11.*closed; %+G11*0.5*g.*G11+G12*0.5*g.*G12+G11.*G12*g.*tanh(beta*phi/2);
end

Fphi=-2/beta*log(2+2*cosh(beta*phi))+phi.^2/g;
Wphi=phi.*exp(-beta*Fphi);
G0=sum(Gphi.*Wphi)/sum(closed.*Wphi);
Sig0=iw0-1./G0;

Z=14+2*cosh(beta);
G0exact=((1/(iw0+1)+1/(iw0-1))*(1+cosh(beta))+12/iw0)/Z;
Sig0exact=iw0-1/G0exact;

end
