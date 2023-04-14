function [DOA]=PAST(X, N_alpha)

[N, T]= size(X);
V_all=[];
W=zeros(N,N_alpha);
% W(1,1)=1;W(3,2)=1;
W(1:N_alpha,1:N_alpha)=eye(N_alpha);
P=eye(N_alpha);
beta=0.97;%%forgetting factor

for tt=1:T
    y1=W'*X(:,tt);
    h=P*y1;
    g=h/(beta+y1'*h);
    P=P-g*h';
    P=(P+P')/2;
    P=P/beta;
    e1=X(:,tt)-W*y1;
    W=W+e1*g';
    V1=W*W';  
    V_all(:,:,tt)=orth(V1);
end

reslu=0.1;
search_area=-90:reslu:90;
for tt=1:T
    if tt==513
       111; 
    end
    Vtt= V_all(:,:,tt);
    Pm(:,tt)=music_linear(Vtt,reslu);
    [DOA(tt,:)]=findmax(search_area,Pm(:,tt),2,0.01);
    %    figure(4);stem(search_area,Pm(:,tt));hold off;
end
