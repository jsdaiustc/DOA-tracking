function theta=findmax(areas,Pm,N,resolve1)
zz=[];
Pm=Pm(:);
d_Pm=diff(Pm);
d_Pm=[-1;d_Pm;1];

Pm_choice=[];

dpos=find(d_Pm>=0);
jj=1;
for ii=1:length(dpos)-1
    if d_Pm(dpos(ii)+1)<=0
      zz(jj)= dpos(ii);
      jj=jj+1;
    end
end

theta=zeros(N,1);

if length(zz)<N
    theta=zeros(N,1);
elseif  length(zz)==N
    theta(1:length(zz))=areas(zz); 
else
    [~,I]=sort(Pm(zz), 'descend');    
    temp_alpha=areas(zz(I)); temp_Pm=Pm(zz(I));

    theta=[];
    theta(1)=temp_alpha(1); Pm_choice(1)=temp_Pm(1);
    ii=2;
    while length(theta)~=N
        if min(abs(theta- temp_alpha(ii)))> resolve1
            theta=[theta; temp_alpha(ii)];
             Pm_choice=[Pm_choice; temp_Pm(ii)];
        end
        ii=ii+1;
    end

end

% [P_max,ind_max]=max(Pm_choice);
% ind_remove=find(Pm_choice< P_max/100);
% theta(ind_remove)= ones(size(ind_remove))* theta(ind_max);



[theta,~]=sort(theta);
theta=theta';


