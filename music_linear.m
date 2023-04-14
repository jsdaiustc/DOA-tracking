function Pm=music_linear(V,reslu)

[M,N_alpha]=size(V);


% [U,T]=eig(R);   
% lam=diag(T);
% [lam,I]=sort(lam,'descend');
% U=U(:,I);    
% U=U(:,N_alpha+1:M);

UUn= eye(M)-V*V';





alfa1=-90:reslu:90;
alfa=alfa1(:)*pi/180.;%????????????????????????????
i=sqrt(-1);
L=(0:M-1)';%(0:M-1)=(0:length(R))
for l=1:length(alfa1)
    
    
%    a=cos(-i*pi*sin(alfa(l))*L);
    
%     if ceil(M/2)- M/2==0  %  M is a even
%        a=[ cos(-pi*sin(alfa(l))*[ (M-1)/2 :-1: 1/2]');
%            sin(-pi*sin(alfa(l))*[1/2:1:(M-1)/2 ]')];
%     else
%         a=[ cos(-pi*sin(alfa(l))*[ (M-1)/2 :-1: 1]');
%             1;
%             sin(-pi*sin(alfa(l))*[1:1:(M-1)/2 ]')];
%     end




a=exp(-i*pi*sin(alfa(l))*L);
tmp=a'*UUn*a;
    
Pm(l)=1/abs(tmp);
    
    
    
    
    
    
end


temp=Pm-min(Pm);

Pm=temp/max(temp);