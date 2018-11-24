Lx=10:2:32;
Ly=10:2:32;

for i=1:length(Lx)
    for j=1:length(Ly)
        S(:,i,j)=soti_qti(Lx(i),Ly(j),0,0);
        %[Lx(i),Ly(j)]
    end
end
dS(1:length(Lx),1:length(Ly))=S(1,:,:)+S(2,:,:)-S(3,:,:);

subplot(2,2,1);
x=[0,2,2,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
%fill(x,y,'g');
hold on;
plot(x,y,'ko','markersize',10,'markerfacecolor','r');
plot([0,2],[1,1],'k');
plot([1,1],[0,2],'k');
text(0.5,1.5,'A');
text(0.5,0.5,'B');
text(1.5,0.5,'C');
text(1.5,1.5,'D');
text(0,2.2,'1');
text(0,-0.2,'2');
text(2,-0.2,'3');
text(2,2.2,'4');
xticks([]);
yticks([]);
xlabel('x');
ylabel('y');
title('(a)');

subplot(2,2,2);
y(1:3,1:length(Ly))=S(1:3,6,:);
plot(Ly,y(1,:),'bo-');
hold on
plot(Ly,y(2,:),'gs-');
plot(Ly,y(3,:),'rv-');
xlim([10,32]);
xlabel('L_y');
ylabel('S_E');
title('(b)');
text(25,27,'S_E(AB)','color','g')
text(24,10,'S_E(AD)','color','b')
text(15,34,'S_E(AC)','color','r')



subplot(2,2,3);
i=3;
plot(Ly,dS(i,:),'b*-');
hold on;
i=6;
plot(Ly,dS(i,:),'g*-');
i=9;
plot(Ly,dS(i,:),'r*-');
i=12;
plot(Ly,dS(i,:),'c*-');
xlim([10,40]);
ylim([0.3,0.8])
xlabel('L_y');
ylabel('\DeltaS_E');
title('(c)');
text(11,0.72,'L_x=14','color','b');
text(19,0.72,'20','color','g');
text(25,0.72,'26','color','r');
text(31,0.72,'32','color','c');

plot([Ly,40],[diag(dS)',dS(length(Lx),length(Ly))],'k--');
plot([Ly,40],[dS(6,:),dS(6,length(Ly))],'k--');

text(32,0.5,'0.2797','color','k')

subplot(2,2,4);
pcolor(Lx,Ly,dS);
shading flat;
colorbar
xticks([10,20,30]);
yticks([10,20,30]);
xlabel('L_x');
ylabel('L_y');
title('(d)');
