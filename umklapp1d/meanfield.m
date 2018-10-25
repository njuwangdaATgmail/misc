function [G0,Sig0,G0exact,Sig0exact]=meanfield(beta)

gap=1;
err=1;
while err>1e-12
    gapnew=tanh(beta*gap/2);
    err=abs(gapnew-gap);
    gap=gapnew;
end

iw0=1i*pi/beta;
G0=0.5*(1/(iw0-gap)+1/(iw0+gap));
Sig0=iw0-1/G0;

Z=14+2*cosh(beta);
G0exact=((1/(iw0+1)+1/(iw0-1))*(1+cosh(beta))+12/iw0)/Z;
Sig0exact=iw0-1/G0exact;

end