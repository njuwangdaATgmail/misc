[~,u,v]=soti(20,20,0,0);
[~,f1]=soti_ky(100,100,1);
E=soti_ky(100,100,0);

subplot(2,2,1);
plot(diag(v),'b.');
xlim([0,1600]);
ylim([-3,3]);
ylabel('energy')
title('(a)');
axes('position',[0.38,0.65,0.1,0.1]);
plot(diag(v),'bo','markersize',4,'markerfacecolor','b');
xlim([798,803]);
ylim([-0.1,0.1]);
hold on
plot([700,900],[0,0],'k:');
xticks([]);
yticks([-0.1,0,0.1]);

subplot(2,2,2);
ldos=zeros(20,20);
for i=799:802
    a=reshape(u(:,i),4,20,20);
    x(1:20,1:20)=sum(abs(a).^2,1);
    ldos=ldos+x;
end
imagesc(ldos);
shading interp;
colorbar;
xticks([]);
yticks([]);
title('(b)');

subplot(2,2,3);
plot((0:100)/100*2,E,'b');
xlim([0,2]);
ylim([-3,3]);
xlabel('k_y');
ylabel('energy');
xticks([0,1,2]);
xticklabels({'0','\pi','2\pi'});
title('(c)');


subplot(2,2,4);
plot((0:100)/100*2,f1,'b');
xlim([0,2]);
ylim([0,1]);
xlabel('k_y');
ylabel('f_E')
xticks([0,1,2]);
xticklabels({'0','\pi','2\pi'});
title('(d)');
axes('position',[0.63,0.23,0.1,0.1]);
x=[0,1,1,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=x+1;
plot(x,y,'k');
axis off