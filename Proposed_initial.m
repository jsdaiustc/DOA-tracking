%Message passing for Group 1
function [DOA,initial]=Proposed_initial(Y,resolution,etc)
[M,N]=size(Y);
search_area_orginal=[-90:resolution:90]';
Phi=exp(-1i*pi*(0:M-1)'*sind(search_area_orginal'))/sqrt(M);
B=-1i*pi*(0:M-1)'*cosd(search_area_orginal').*Phi;
L=length(search_area_orginal);
search_area=repmat(search_area_orginal,1,N);
Phi=repmat(Phi,1,1,N);
B=repmat(B,1,1,N);

maxiter=200;
tol=1e-5;
converged = false;
iter = 0;
alpha_inv=1;
gamma_inv=ones(L,1);

% GAMP initialization
mu=zeros(L,N);
V_x=ones(L,N);
hold_max=10^10;
hold_min=10^(-10);
mu_y  =  zeros(M,N);
A2 =  1/M;
rho=0.7;     % parameter of damping

while ~converged
    %% calculate mu and Sigma
    alpha=1/alpha_inv;
    Temp_ay=zeros(L,N);
    mu_tt_old=mu;
    V_x_old=V_x;
    for nn=1:N
        h(:,nn)= Phi(:,:,nn)*mu(:,nn);
    end
    V_p     = (ones(M,1)* A2)* sum(V_x,1);         %EQ.49
    V_p     = min(max(V_p,hold_min),hold_max);
    mu_p    = h-  V_p.*mu_y;                       %EQ.50
    V_h     = (V_p*alpha)./(V_p+alpha);            %EQ.52
    mu_h    = V_h .*(Y./alpha+mu_p./V_p);          %EQ.53
    mu_y    = (mu_h -mu_p)./V_p;                   %EQ.48
    V_y     = (1-V_h./V_p)./V_p;                   %EQ.47
    V_y     = min(max(V_y,hold_min),hold_max);
    V_acute_inv = ones(L,1)* (A2*sum(V_y,1));
    V_acute = 1./V_acute_inv;                      %EQ.43
    for nn=1:N
        Temp_ay(:,nn)= Phi(:,:,nn)'* mu_y(:,nn);
    end
    s_acute = mu + V_acute.*Temp_ay;                         %EQ.44
    mu = (s_acute.*V_acute_inv)./(gamma_inv + V_acute_inv);  %EQ.46
    V_x = 1./(gamma_inv+V_acute_inv);                        %EQ.45
    mu= (1-rho)*mu + (rho)*mu_tt_old;                        %damping
    V_x= (1-rho)*V_x + rho*V_x_old;                          %damping

    
    %% update alpha0
    alpha_rho=0.99;
    alpha_inv_old=alpha_inv;
    resid_all =  Y- mu_h;
    alpha_inv=(N*M)/( norm(resid_all, 'fro')^2 +   sum(V_h(:)) );   %EQ.61
    alpha_inv=alpha_rho*alpha_inv_old + (1-alpha_rho)*alpha_inv;    %damping
    
    %% update gamma
    gamma_inv_old=gamma_inv;
    mu2= mu.*conj(mu);
    tc = sum( mu2 +  V_x,2 );
    gamma_inv=N./tc;                                                %EQ.64
    
    %% update grid
    search_area_final=search_area;
    Pm=sum(mu2,2);
    [~,sort_ind]=sort(Pm, 'descend');
    active_set=sort_ind(1:etc);
    P=zeros(etc,etc,N); r=zeros(etc,N);
    for nn=1:N
        BHB = B(:,active_set,nn)' * B(:,active_set,nn);
        P(:,:,nn) = real( conj(BHB) .* ( (mu(active_set,nn) * mu(active_set,nn)') + diag(V_x(active_set,nn))  )  );
        Tes=diag(V_x(:,nn));
        v2 = real(diag(B(:,active_set,nn)' * Phi(:,:,nn) *Tes(:,active_set) ));
        r(:,nn) =  real(conj(mu(active_set,nn)) .* (B(:,active_set,nn)' * (Y(:,nn) - Phi(:,:,nn)  * mu(:,nn)) )   )  - v2;
    end
    temp_grid1=sum(r,2)./(diag(sum(P,3))+ 1e-20);
    temp_grid1=temp_grid1*180/pi;
    theld=resolution/20*0.95^(iter);
    ind_small=find(abs(temp_grid1)<theld);
    temp_grid1(ind_small)=sign(temp_grid1(ind_small))*theld;
    ind_unchang=find (abs(temp_grid1)>resolution);
    temp_grid1(ind_unchang)=sign(temp_grid1(ind_unchang)) * resolution/20;
    for nn=1:N
        search_area(active_set,nn)=search_area(active_set,nn) + temp_grid1;
    end
    P_d= P(:,:,1) + P(:,:,3);
    v_d= r(:,3)-r(:,1);
    temp_d= v_d./(diag(P_d) + 1e-20);
    ind_unchang=find (abs(temp_d)>0.1);  temp_d(ind_unchang)=0.1;
    temp_d=temp_d/max(abs(temp_d))*0.01*0.995^(iter);
    weight=[-1,0,1];
    for nn=1:N
        search_area(active_set,nn)=search_area(active_set,nn) + weight(nn)*temp_d;
    end
    for nn=1:N
        Phi(:,active_set,nn)=exp(-1i*pi*(0:M-1)'*sind(search_area(active_set,nn)'))/sqrt(M);
        B(:,active_set,nn)=-1i*pi*(0:M-1)'*cosd(search_area(active_set,nn)').*Phi(:,active_set,nn);
    end
    
    %% check the convergence
    erro=norm(gamma_inv_old - gamma_inv)/norm(gamma_inv_old);
    if  iter >= maxiter || erro<tol
        converged = true;
    end
    iter = iter + 1;
end

dis=0.1;
Pm=sum(mu.*conj(mu),2);
for nn=1:N
    search_ee=search_area_final(:,nn);
    [search_ee,ind_ee]=sort(search_ee,'ascend');
    Pm_ee=Pm(ind_ee);
    insert=(search_ee(1:end-1)+search_ee(2:end))/2;
    search_area_2=zeros(length(search_ee)*2-1,1);
    search_area_2(1:2:end)=search_ee;
    search_area_2(2:2:end)=insert;
    Pm_2=zeros(length(Pm_ee)*2-1,1);
    Pm_2(1:2:end)=Pm_ee;
    [DOA(nn,:)]=findmax(search_area_2,Pm_2,2,dis);
end

initial.gamma_inv=gamma_inv(:,end);
initial.search_area=search_area;
initial.Phi=Phi;
initial.B=B;