function[h]=hk_(kx,ky)
% absolute coordinate is used

t1s=1.79;
t1p=-0.45*t1s;
t2s=0.04*t1s;
t2p=-0.15*t1s;
soc=0.35*t1s;

theta1=(0:2)*2*pi/3;  % nearest neighbour bond angles from a to b
theta2=(1:6)*pi/3-pi/6;  % next nearest neighbour bond angles

eikr1=exp(1i*kx*cos(theta1)+1i*ky*sin(theta1));
eikr2=exp(1i*kx*sqrt(3)*cos(theta2)+1i*ky*sqrt(3)*sin(theta2));

h=zeros(4,4); % basis defined as [xa,ya,xb,yb];

% nearest neighbour hopping from a to b
h(1,3) = t1s*sum( cos(theta1).*cos(theta1).*eikr1 ) + t1p*sum( sin(theta1).*sin(theta1).*eikr1 );
h(2,4) = t1s*sum( sin(theta1).*sin(theta1).*eikr1 ) + t1p*sum( cos(theta1).*cos(theta1).*eikr1 );
h(1,4) = t1s*sum( cos(theta1).*sin(theta1).*eikr1 ) - t1p*sum( sin(theta1).*cos(theta1).*eikr1 );
h(2,3) = t1s*sum( sin(theta1).*cos(theta1).*eikr1 ) - t1p*sum( cos(theta1).*sin(theta1).*eikr1 );
h(3:4,1:2)=h(1:2,3:4)';

% next nerest neighbour hopping from a to a
h(1,1) = t2s*sum( cos(theta2).*cos(theta2).*eikr2 ) + t2p*sum( sin(theta2).*sin(theta2).*eikr2 );
h(2,2) = t2s*sum( sin(theta2).*sin(theta2).*eikr2 ) + t2p*sum( cos(theta2).*cos(theta2).*eikr2 );
h(1,2) = t2s*sum( cos(theta2).*sin(theta2).*eikr2 ) - t2p*sum( sin(theta2).*cos(theta2).*eikr2 );
h(2,1) = t2s*sum( sin(theta2).*cos(theta2).*eikr2 ) - t2p*sum( cos(theta2).*sin(theta2).*eikr2 );

% next nearest neighbour hopping from b to b
h(3:4,3:4)=h(1:2,1:2);  % check it!

% soc
h(1,2)=h(1,2)-soc*(-1i);
h(2,1)=h(2,1)-soc*(1i);
h(3,4)=h(3,4)-soc*(-1i);
h(4,3)=h(4,3)-soc*(1i);

h=(h+h')/2;

end
