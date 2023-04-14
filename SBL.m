function DOA_full=SBL(Y,resolution,etc,ex)
[M,T]=size(Y);
Y_norm=norm(Y,'fro')/sqrt(M*T);
Y=Y/Y_norm;

T=T/ex;
search_area_orginal=[-90:resolution:90]';
search_area=search_area_orginal*ones(1,T);

for t=1:T
    Phi(:,:,t)=exp(-1i*pi*(0:M-1)'*sind(search_area(:,t)'))/sqrt(M);
    B(:,:,t)=-1i*pi*(0:M-1)'*cosd(search_area(:,t)').*Phi(:,:,t);
end

[K_hat,~]=size(search_area);

a=1e-20;
b=a;
maxiter=200;
tol=1e-5;
converged = false;
iter = 0;
alpha0=ones(T,1);
D=ones(K_hat,T);
gamma=ones(K_hat,T)*1;

mu=zeros(K_hat, T*ex);
Sigma=ones(K_hat,K_hat, T);

while ~converged
    %%%%%%%%%%%%%%%%%%%%%calculate mu and Sigma
    

    mu_old=mu;
    Sigma_old=Sigma;
    for tt=1:T
        D(:,tt)=gamma(:,tt) ;
        diag_inv=diag(  1./D(:,tt)  );
        Phi_diag(:,:,tt) = Phi(:,:,tt) *  (diag_inv);
        V_temp(:,:,tt)= 1/alpha0(tt)*eye(M) + Phi_diag(:,:,tt) * Phi(:,:,tt)';
        Sigma(:,:,tt) = (diag_inv) -Phi_diag(:,:,tt)' * (V_temp(:,:,tt)\Phi_diag(:,:,tt));%woodbury
        mu(:,(ex*(tt-1)+1):ex*tt)=  alpha0(tt)*Sigma(:,:,tt)*Phi(:,:,tt)'*Y(:,(ex*(tt-1)+1):ex*tt);
    end
%     mu=rho1*mu_old + (1-rho1)*mu;
%     for tt=1:T
%        Sigma(:,:,tt)=rho1*Sigma_old(:,:,tt) + (1-rho1)*Sigma(:,:,tt);
%     end
    
    
    
    %% update alpha0
%     alpha0_old=alpha0;
%      for tt=1:T
%         resid=Y(:,(ex*(tt-1)+1):ex*tt) - Phi(:,:,tt)*mu(:,(ex*(tt-1)+1):ex*tt);
%         be_b= norm(resid, 'fro')^2+ ex* real( trace(Phi(:,:,tt)*Sigma(:,:,tt)*Phi(:,:,tt)'));
%         alpha0(tt)=(a+ex*M)/(b+ be_b);
%     end
%     alpha_rho=0.98;
%     alpha0=alpha_rho*alpha0_old + (1-alpha_rho)*alpha0;
%     
    
     alpha_rho=0.98;
     alpha0_old=alpha0;
     for tt=1:T
        resid=Y(:,(ex*(tt-1)+1):ex*tt) - Phi(:,:,tt)*mu(:,(ex*(tt-1)+1):ex*tt);
        be_b(tt)= norm(resid, 'fro')^2+ ex* real( trace(Phi(:,:,tt)*Sigma(:,:,tt)*Phi(:,:,tt)'));
%         alpha0(tt)=(a+ex*M)/(b+ be_b);
     end
    alpha0= (a+ex*M*T)/(b+ sum(be_b));
    alpha0=alpha0*ones(T,1);
    alpha0=alpha_rho*alpha0_old + (1-alpha_rho)*alpha0   ;
    
    
    
    
    
    
    
    %% update gamma
    mu2= mu.*conj(mu);
    temp=zeros(K_hat,T);
    for tt=1:T
        for i=1:ex
            temp(:,tt)=temp(:,tt)+  mu2(:,(ex*(tt-1)+i)) +  real(  diag( Sigma(:,:,tt) )  );
        end
    end
    
    tc=temp;

    c_k=a+ex;
    d_k=b+(tc);
    gamma=c_k./d_k;

    
    
 
    
    %% update grid
    
    %     figure(5);stem(search_area(:,end/2),  1./gamma(:,end/2))
    search_area_old=search_area;
    active_set=zeros(etc,T);
    
    for tt=1:T
         Pm(:,tt)=sum(mu(:,(ex*(tt-1)+1):ex*tt).*conj(mu(:,(ex*(tt-1)+1):ex*tt)),2);
        [~,sort_ind]=sort(Pm(:,tt), 'descend');
        idx=sort_ind(1:etc);
        active_set(:,tt)=idx;
        
        BHB = B(:,idx,tt)' * B(:,idx,tt);
        P = real( conj(BHB) .* ((mu(idx,(ex*(tt-1)+1):ex*tt) * mu(idx,(ex*(tt-1)+1):ex*tt)') + ex*  Sigma(idx,idx,tt)   )  );
        v2= ex*  real(diag(B(:,idx,tt)' * Phi(:,:,tt) * Sigma(:,idx,tt)));
        v = sum(   real(conj(mu(idx,(ex*(tt-1)+1):ex*tt)) .* (B(:,idx,tt)' * (Y(:,(ex*(tt-1)+1):ex*tt) - Phi(:,:,tt)  * mu(:,(ex*(tt-1)+1):ex*tt)) )   ),2) -   v2;
        
        temp_grid1=v./(diag(P)+ 1e-20);
%         temp_grid1=pinv(P)*v;
        temp_grid1=temp_grid1*180/pi;
        theld=resolution/20*0.95^(iter);
        ind_small=find(abs(temp_grid1)<theld);
        temp_grid1(ind_small)=sign(temp_grid1(ind_small))*theld;
        ind_unchang=find (abs(temp_grid1)>resolution);
        temp_grid1(ind_unchang)=sign(temp_grid1(ind_unchang)) * resolution/20;
        search_area(idx,tt)=search_area(idx,tt) + temp_grid1;
     end


    
    
    for tt=1:T
        Phi(:,active_set(:,tt),tt)=exp(-1i*pi*(0:M-1)'*sind(search_area(active_set(:,tt),tt)'))/sqrt(M);
        B(:,active_set(:,tt),tt)=-1i*pi*(0:M-1)'*cosd(search_area(active_set(:,tt),tt)').*Phi(:,active_set(:,tt),tt);
    end
    
 
    
    %     erro=norm(D - diag_last)/norm(diag_last);
    if  iter >= maxiter
        converged = true;
    end
    
    
    iter = iter + 1;
end

for t=1:T
    Pm=sum(mu(:,(ex*(t-1)+1):ex*t).*conj(mu(:,(ex*(t-1)+1):ex*t)),2);
    %     plot(search_area(:,t),Pm)
    [DOA(t,:)]=findmax(search_area_old(:,t),Pm,2,resolution);
end
DOA_full=DOA;
% DOA_full=zeros( T*ex, size(DOA,2));
% for ii=1:size(DOA,2)
%     xi=1:1:T*ex;
%     DOA_full(:,ii)=interp1([1:ex:T*ex],DOA(:,ii),xi, 'spline');
% end
