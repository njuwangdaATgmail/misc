[~,~,~,~,~,~,~,f4,f5]=soti(21,21,1,1);
[~,~,~,f1,f2,~,f3]=soti(20,20,1,1);

subplot(2,2,1);
plot(f1,'b.');
xlim([0,800]);
ylim([0,1]);
xticks([0,400,800]);
ylabel('f_E');
title('(a)');
axes('position',[0.17,0.75,0.1,0.1]);
x=[0,1,1,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=x+1;
plot(x,y,'k');
axis off


subplot(2,2,2);
plot(f3,'b.');
xlim([0,400]);
ylim([0,1]);
ylabel('f_E');
xticks([0,200,400]);
title('(b)');
axes('position',[0.82,0.65,0.1,0.1]);
plot(f3,'bo','markersize',3,'markerfacecolor','b');
hold on;
plot([100,300],[0.5,0.5],'k:')
xlim([194,207]);
ylim([0.2,0.8]);
xticks([]);
yticks([0.2,0.5,0.8]);
axes('position',[0.6,0.75,0.1,0.1]);
x=[0,1,1,0,0];
y=[0,0,1,1,0];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=x*2;
y=y*2;
plot(x,y,'k');
axis off

subplot(2,2,3);
plot(f4,'b.');
xlim([0,900]);
ylim([0,1]);
xticks([0,400,800]);
ylabel('f_E');
title('(c)');
axes('position',[0.38,0.17,0.1,0.1]);
plot(f4,'bo','markersize',3,'markerfacecolor','b');
hold on;
plot([400,500],[0.5,0.5],'k:')
xlim([438,447]);
ylim([0.45,0.55]);
xticks([]);
yticks([0.45,0.5,0.55]);
axes('position',[0.17,0.27,0.1,0.1]);
x=[0,1,2,1,0];
y=[1,0,1,2,1];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=[0,2,2,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
axis off

subplot(2,2,4);
plot(f5,'b.');
xlim([0,450]);
ylim([0,1]);
xticks([0,200,400]);
ylabel('f_E');
title('(d)');
axes('position',[0.6,0.27,0.1,0.1]);
x=[0.5,1,1.5,1,0.5];
y=[1,0,1,2,1];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=[0,2,2,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
axis off
