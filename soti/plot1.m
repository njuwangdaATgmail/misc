[S,u,v,f1,f2,f3,f4,f5]=soti(21,21,1,1);
[S,u,v,f1,f2,f3]=soti(20,20,1,1);

subplot(2,2,1);
plot(f1,'.');
xlim([0,800]);
ylim([0,1]);
title('(a)');
axes('position',[0.15,0.75,0.1,0.1]);
x=[0,1,1,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=x+1;
plot(x,y,'k');
axis off


subplot(2,2,2);
plot(f3,'.');
xlim([0,400]);
ylim([0,1]);
title('(b)');
axes('position',[0.82,0.65,0.1,0.1]);
plot(f3,'.');
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
plot(f4,'.');
xlim([0,900]);
ylim([0,1]);
title('(c)');
axes('position',[0.38,0.2,0.1,0.1]);
plot(f4,'.');
hold on;
plot([400,500],[0.5,0.5],'k:')
xlim([438,447]);
ylim([0.45,0.55]);
xticks([]);
yticks([0.45,0.5,0.55]);
axes('position',[0.15,0.25,0.1,0.1]);
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
plot(f5,'.');
xlim([0,450]);
ylim([0,1]);
title('(d)');
axes('position',[0.6,0.25,0.1,0.1]);
x=[0.5,1,1.5,1,0.5];
y=[1,0,1,2,1];
plot(x,y,'k');
fill(x,y,'g');
hold on
x=[0,2,2,0,0];
y=[0,0,2,2,0];
plot(x,y,'k');
axis off
