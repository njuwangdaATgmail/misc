function[Ecut,Ecutbdg,kcut]=band(mu,gap)

uca=[3/2,sqrt(3)/2];
ucb=[0,sqrt(3)];

ucvol=3*sqrt(3)/2;

ucka=2*pi*sqrt(3)/ucvol*[1,0];
uckb=2*pi*sqrt(3)/ucvol*[-1/2,sqrt(3)/2];

K=ucka/3+2*uckb/3;
Kp=2*ucka/3+uckb/3;
M=ucka/2+uckb/2;

Nk=500;
kcut=zeros(3*Nk+1,2);
kcut(1:Nk,1)=(0:Nk-1)/Nk*M(1);
kcut(1:Nk,2)=(0:Nk-1)/Nk*M(2);
kcut(Nk+1:2*Nk,1)=M(1)+(0:Nk-1)/Nk*(K(1)-M(1));
kcut(Nk+1:2*Nk,2)=M(2)+(0:Nk-1)/Nk*(K(2)-M(2));
kcut(2*Nk+1:3*Nk,1)=K(1)+(0:Nk-1)/Nk*(-K(1));
kcut(2*Nk+1:3*Nk,2)=K(2)+(0:Nk-1)/Nk*(-K(2));
kcut(3*Nk+1,1:2)=0;

Ecut=zeros(4,3*Nk+1);
Ecutbdg=zeros(8,3*Nk+1);
for i=1:length(kcut)
  h=hk(kcut(i,1),kcut(i,2));
  Ecut(1:4,i)=eig(h);
  H=bdgk(kcut(i,1),kcut(i,2),mu,gap);
  Ecutbdg(1:8,i)=eig(H);
end

end
