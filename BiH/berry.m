function[lam,Ecut]=berry(b,mu,gap)

ucvol=3*sqrt(3)/2;

ucka=2*pi*sqrt(3)/ucvol*[1,0];
uckb=2*pi*sqrt(3)/ucvol*[-1/2,sqrt(3)/2];

Nk=20000;
ka=zeros(Nk,2);
ka(:,1)=(1:Nk)/Nk*ucka(1);
ka(:,2)=(1:Nk)/Nk*ucka(2);
kb=b*uckb;

P1=eye(8);
P2=eye(8);
P3=eye(8);
P4=eye(8);
P5=eye(8);
P6=eye(8);
P7=eye(8);
P8=eye(8);

Ecut=zeros(8,Nk);
for i=1:Nk
    Hk=bdgk(ka(i,1)+kb(1),ka(i,2)+kb(2),mu,gap);
    [u,~]=eig(Hk);
    u(3:4,:)=u(3:4,:)*exp(1i*(ka(i,1)+kb(1)));
    u(7:8,:)=u(7:8,:)*exp(1i*(ka(i,1)+kb(1)));
    Ecut(:,i)=eig(Hk);
    P1=P1*u(:,1)*u(:,1)';
    P2=P2*u(:,2)*u(:,2)';
    P3=P3*u(:,3)*u(:,3)';
    P4=P4*u(:,4)*u(:,4)';
    P5=P5*u(:,5)*u(:,5)';
    P6=P6*u(:,6)*u(:,6)';
    P7=P7*u(:,7)*u(:,7)';
    P8=P8*u(:,8)*u(:,8)';
end


lam=angle([trace(P1),trace(P2),trace(P3),trace(P4),trace(P5),trace(P6),trace(P7),trace(P8)]);

end