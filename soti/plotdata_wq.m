% % fig3
% load soti;
% subplot(2,1,1);
% plot(A,oS,'b:')
% hold on
% plot(A,oS,'*m')
% hold off
% text(20,0.75,'\DeltaS=S(AB)+S(AD)-S(AC)','fontsize',15)
% text(33,0.64,'Universal','color','red','fontsize',15)
% text(33,0.54,'value','color','red','fontsize',15)
% text(33,0.44,'0.2797','color','red','fontsize',15)
% % annotattext('arrow',[0.5,0.37],[0.5,0.64])
% ylim([0.3,0.8])
% hold off
% xlabel('L_y')
% ylabel('\DeltaS')
% title('(a)','fontsize',15);
% 
% subplot(2,1,2);
% pcolor(10:2:32,10:2:32,oS)
% shading interp
% colorbar
% xlabel('L_x')
% ylabel('L_y')
% title('(b)','fontsize',15);



% fig2
m1=[0,2,2,0,0];m2=[0,0,2,2,0];m3=[0,2];m4=[1,1];
subplot(2,1,1);
plot(m1,m2)
fill(m1,m2,'g')
axis off
hold on
plot(m3,m4,'r',m4,m3,'m')
plot(m1,m2,'ko','markersize',10,'markerfacecolor','k')
hold off
title('(a)','fontsize',15)
text(0.5,1.5,'A','fontsize',20)
text(0.5,0.5,'B','fontsize',20)
text(1.5,0.5,'C','fontsize',20)
text(1.5,1.5,'D','fontsize',20)
text(1,-0.2,'L_x')
text(-0.1,1.7,'L_y')

load soti;
subplot(2,1,2);
plot(A,S1(6,:),'b:',A,S2(6,:),'b:',A,S3(6,:),'b:')
hold on
plot(A,S1(6,:),'kx',A,S2(6,:),'ro',A,S3(6,:),'g+')
hold off
xlabel('L_y')
ylabel('Entanglement entropy')
text(20,37,'S(AC)','fontsize',15)
text(30,30,'S(AB)','fontsize',15)
text(27,10,'S(AD)','fontsize',15)
title('(b)','fontsize',15)



% fig1
figure
load soti
subplot(2,2,1);
plot(E,'.')
set(gca,'xtick',[0:400:1600],'ytick',[-3:1:3])
ylabel('Entanglement entropy')
title('(a)','fontsize',15)
axes('position',[0.35,0.8,0.15,0.1])
ylim([-0.01,0.01])
plot(E(799:802),'.')
ylim([-0.01,0.01])
xlim([0,5])
set(gca,'xtick',[])

subplot(2,2,2);
imagesc(a3);
shading interp
colorbar
xlabel('L_x')
ylabel('L_y')
title('(b)','fontsize',15)

subplot(2,2,3);
plot(ky_array,E1,'b.')
text(2.8,-3.5,'k_y')
ylabel('Energy')
set(gca,'xtick',[],'ytick',-3:1:3)
xlim([0,2*pi])
text(pi,-3.2,'\pi','fontsize',15)
text(2*pi,-3.2,'2\pi','fontsize',15)
title('(c)','fontsize',15)

subplot(2,2,4);
plot(ky_array,lam,'.')
text(2.8,-.2,'k_y')
ylabel('Entanglement entropy')
set(gca,'xtick',[])
xlim([0,2*pi])
ylim([-0.1,1.1])
text(pi,-.12,'\pi','fontsize',15)
text(2*pi,-.12,'2\pi','fontsize',15)
title('(d)','fontsize',15)
% % 
