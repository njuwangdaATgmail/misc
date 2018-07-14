%[S,u,v,f1,f2,f3]=soti(20,20,0,0);

subplot(2,2,1);
plot(diag(v),'.');
xlim([0,1600]);
title('(a)');

subplot(2,2,2);
ldos=zeros(20,20);
for i=799:802
    a=reshape(u(:,i),4,20,20);
    x(1:20,1:20)=sum(abs(a).^2,1);
    ldos=ldos+x;
end
imagesc(ldos);
shading interp;
colorbar

subplot(2,2,3);
