function[Hk]=bdgk(kx,ky,mu,gap)

Hk=zeros(8);

h=hk(kx,ky);
Hk(1:4,1:4)=h-mu*eye(4);
h=hk(-kx,-ky);
Hk(5:8,5:8)=-conj(h-mu*eye(4));

Hk(1,6)=-1i*gap;
Hk(2,5)=1i*gap;
Hk(3,8)=1i*gap;
Hk(4,7)=-1i*gap;

Hk(5:8,1:4)=Hk(1:4,5:8)';

Hk=(Hk+Hk')/2;

end
