function[G0,Sig0,G0exact,Sig0exact]=uniformphi(beta)

dphi=0.01;
phi=0:dphi:100*sqrt(1/beta);  % |phi|^2>>T*g, where g=1
iw0=1i*pi/beta;
Gphi=0.5*(1./(iw0-phi)+1./(iw0+phi));
Fphi=-2/beta*log(2+2*cosh(beta*phi))+phi.^2;
Wphi=phi.*exp(-beta*Fphi);
G0=sum(Gphi.*Wphi)/sum(Wphi);
Sig0=iw0-1./G0;

Z=14+2*cosh(beta);
G0exact=((1/(iw0+1)+1/(iw0-1))*(1+cosh(beta))+12/iw0)/Z;
Sig0exact=iw0-1/G0exact;

end