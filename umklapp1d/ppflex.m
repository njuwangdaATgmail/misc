function[G,self,chi,T11,wn_array,vn_array,k_array]=ppflex(g3,g4,T)

mixing=1;
Nw=1000;
Nk=100;

% wn_array(1<=i<=2*Nw)=(2*i-2*Nw-1)*pi*T;
% vn_array(1<=i<=4*Nw-1)=(2*i-4*Nw)*pi*T;

wn_array=(-2*Nw+1:2:2*Nw-1)*pi*T;
vn_array=(-4*Nw+2:2:4*Nw-2)*pi*T;
k_array=(0:Nk-1)*2*pi/Nk;

ek=-2*cos(k_array);

G=zeros(length(wn_array),length(k_array));
self=G;
selfnew=G;
chi=zeros(length(vn_array),length(k_array));
T11=chi;

error=1;
iter=0;
while error>1e-6
    
    iter=iter+1;

    % get G
    for iw=1:length(wn_array)
        for ik=1:length(k_array)
            G(iw,ik)=1/(1i*wn_array(iw)-ek(ik)-self(iw,ik));
        end
    end

    % get chi
    for iv=1:length(vn_array)
        for iq=1:length(k_array)
            if iv<2*Nw
                seqw=1:iv;
                seqw2=iv:-1:1;
            else
                seqw=1+iv-2*Nw:2*Nw;
                seqw2=2*Nw:-1:1+iv-2*Nw;
            end
            seqk=1:Nk;
            seqk2=mod((Nk+1:-1:2)+iq-2,Nk)+1;
            chi(iv,iq)=sum(sum(G(seqw,seqk).*G(seqw2,seqk2)))*T/length(k_array);
        end
    end
    
    % get T-matrix
    T0v=inv([-g4,-g3;-g3,-g4]);
    for iv=1:length(vn_array)
        for iq=1:length(k_array)
            Tmat=inv(T0v-chi(iv,iq)*eye(2));
            T11(iv,iq)=Tmat(1,1);
        end
    end

    % get self energy
    for iw=1:length(wn_array)
        for ik=1:length(k_array)
            seqv=iw:2*Nw+iw-1;
            seqw2=1:2*Nw;
            seqq=mod((ik:Nk+ik-1)-1,Nk)+1;
            seqk2=1:Nk;
            selfnew(iw,ik)=-sum(sum(T11(seqv,seqq).*G(seqw2,seqk2)))*T/length(k_array);
        end
    end

    error=max(max(abs(selfnew-self)));
    disp([iter,error]);
    self=self+mixing*(selfnew-self);

end

end