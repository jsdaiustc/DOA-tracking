function DOA=BSPLR(Y,resolution,etc,N_alpha)
[M,T]=size(Y);
norm_y=norm(Y,'fro')/sqrt(M*T);
Y=Y/norm_y;

search_area=[-90:resolution:90]';
pos_all=0:M-1;
A_fixed=exp(-1i*pi*pos_all'*sin(search_area'*pi/180) )/sqrt(M);
B_fixed=-1i*pi*pos_all'*cos( search_area'*pi/180 ).*A_fixed;

K_hat=length(search_area);
A_vary=repmat(A_fixed,1,1,T);
B_vary=repmat(B_fixed,1,1,T);
search_area2=repmat(search_area,1,T);





%% initialization
converged = false;
iter = 0;
maxiter=200;
tol=1e-5;
alpha0=1;
delta=ones(K_hat,1);
beta=ones(K_hat,T);

%% calculate mu and Sigma
while ~converged
    for t=1:T
        Phi_tt=[A_fixed,A_vary(:,:,t)];
        Xi=[delta;beta(:,t)];
        Phi_delta = Phi_tt *  diag(Xi);
        V_temp=alpha0*eye(M) + Phi_delta * Phi_tt';
        Sigma(:,:,t) = diag(Xi) -Phi_delta' * (V_temp \Phi_delta);
        %         mu(:,t) =Phi_delta'*(V_temp \Y(:,t));
        mu(:,t) = (1/alpha0) * (  Sigma(:,:,t) * (Phi_tt' * Y(:,t))  );
    end
    
    
    
    %% update delta
    delta_last = delta;
    sum_mu=sum(mu(1:K_hat,:).*conj(mu(1:K_hat,:)),2);
    temp=sum_mu + real(diag( sum(Sigma(1:K_hat,1:K_hat,:),3) ));
    delta=    ( real(temp)  )/T;
    
    %% update beta
    for t=1:T
        beta(:,t)=mu(K_hat+1:end,t).*conj(mu(K_hat+1:end,t))    +     real(diag(Sigma(K_hat+1:end,K_hat+1:end,t)));
    end
    
    %% update alpha
    temp2=0;
    for t=1:T
        Phi_tt=[A_fixed,A_vary(:,:,t)];
        temp2= temp2+real(trace(Phi_tt*Sigma(:,:,t)*Phi_tt'));
        resid(:,t)=Y(:,t)-Phi_tt*mu(:,t);
    end
    alpha_rho=0.99;
    alpha0_old=alpha0;
    alpha0=( norm(resid, 'fro')^2  +   temp2   )/( M*T );
    alpha0=1./(      alpha_rho*(1./alpha0_old) + (1-alpha_rho)*(1./alpha0)     );
    
    
    
    %% grid refine 1
    mu_fixed=mu(1:K_hat,:);
    Sigma_fixed=Sigma(1:K_hat,1:K_hat,:);
    mu_vary=mu(1+K_hat:end,:);
    Sigma_vary=Sigma(1+K_hat:end,1+K_hat:end,:);
    Pm_fixed=sum(mu_fixed.*conj(mu_fixed),2);
    [~,sort_ind_fixed]=sort(Pm_fixed, 'descend');    %降序排列
    idx_fixed=sort_ind_fixed(1:etc);
    BHB_fixed = B_fixed(:,idx_fixed)' * B_fixed(:,idx_fixed);
    P_temp1= sum(Sigma_fixed(idx_fixed,idx_fixed,:),3);
    P_fixed = real( conj(BHB_fixed) .* ((mu_fixed(idx_fixed,:) * mu_fixed(idx_fixed,:)') +   P_temp1   )  );
    v_temp1=real(diag(B_fixed(:,idx_fixed)' * A_fixed * sum(Sigma_fixed(:,idx_fixed,:),3)    ));  %    real(diag(B(:,idx)' * A * Tes(:,idx)   ));
    for t=1:T
        Y_fix(:,t)=Y(:,t)-A_vary(:,:,t) * mu_vary(:,t);
    end
    v_fixed = sum( real(conj(mu_fixed(idx_fixed,:)) .* (B_fixed(:,idx_fixed)' * (Y_fix - A_fixed * mu_fixed ) ) ),2) -   v_temp1;
    temp_grid_fixed=v_fixed./diag(P_fixed);
    temp_grid_fixed=temp_grid_fixed'*180/pi;
    theld=resolution/20*0.95^(iter);
    ind_small=find(abs(temp_grid_fixed)<theld);
    temp_grid_fixed(ind_small)=sign(temp_grid_fixed(ind_small))*theld;
    ind_unchang=find (abs(temp_grid_fixed)>resolution);
    temp_grid_fixed(ind_unchang)=sign(temp_grid_fixed(ind_unchang)) * resolution/20;
    search_area(idx_fixed)=search_area(idx_fixed) + temp_grid_fixed';
    A_fixed(:,idx_fixed)=exp(-1i*pi*pos_all'*sin(search_area(idx_fixed)'*pi/180))/sqrt(M);
    B_fixed(:,idx_fixed)=-1i*pi*pos_all'*cos(search_area(idx_fixed)'*pi/180).*A_fixed(:,idx_fixed);
    
    if iter>20
        %% grid refine 2
        for t=1:T
            Pm_vary=sum(mu_vary(:,t).*conj(mu_vary(:,t)),2);
            [~,sort_ind_vary]=sort(Pm_vary, 'descend');    %降序排列
            idx_vary=sort_ind_vary(1:etc);
            BHB_vary = B_vary(:,idx_vary,t)' * B_vary(:,idx_vary,t);
            P_temp2=  Sigma_vary(idx_vary,idx_vary,t);
            P_vary = real( conj(BHB_vary) .* ((mu_vary(idx_vary,t) * mu_vary(idx_vary,t)') +   P_temp2   )  );
            v_temp2= real(diag(B_vary(:,idx_vary,t)' * A_vary(:,:,t) * Sigma_vary(:,idx_vary,t)));  %    real(diag(B(:,idx)' * A * Tes(:,idx)   ));
            v_vary = sum( real(conj(mu_vary(idx_vary,t)) .* (B_vary(:,idx_vary,t)' * (Y(:,t)-A_fixed*mu_fixed(:,t) - A_vary(:,:,t) * mu_vary(:,t) ) ) ),2) -   v_temp2;
            temp_grid_vary=v_vary./diag(P_vary);
            temp_grid_vary=temp_grid_vary*180/pi;
            theld=resolution/20*0.95^(iter);
            ind_small=find(abs(temp_grid_vary)<theld);
            temp_grid_vary(ind_small)=sign(temp_grid_vary(ind_small))*theld;
            ind_unchang=find (abs(temp_grid_vary)>resolution);
            temp_grid_vary(ind_unchang)=sign(temp_grid_vary(ind_unchang)) * resolution/20;
            search_area2(idx_vary,t)=search_area2(idx_vary,t) + temp_grid_vary;
            A_vary(:,idx_vary,t)=exp(-1i*pi*pos_all'*sin(search_area2(idx_vary,t)'*pi/180))/sqrt(M);
            B_vary(:,idx_vary,t)=-1i*pi*pos_all'*cos(search_area2(idx_vary,t)'*pi/180).*A_vary(:,idx_vary,t);
        end
    end
    
    
    %% stopping criteria
    erro=norm(delta - delta_last)/norm(delta_last);
    if erro < tol || iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
    
end
dis=1;
mu_fixed=mu(1:K_hat,:);
mu_vary=mu(1+K_hat:end,:);
[search_area,ind_sort]=sort(search_area,'ascend');
mu_fixed=mu_fixed(ind_sort,:);
insert=(search_area(1:end-1)+search_area(2:end))/2;
search_area_sort=zeros(length(search_area)*2-1,1);
search_area_sort(1:2:end)=search_area;
search_area_sort(2:2:end)=insert;
mu_fixed_2=zeros(size(mu_fixed,1)*2-1,T);
mu_fixed_2(1:2:end,:)=mu_fixed;
Pm_fixed=mean(mu_fixed_2.*conj(mu_fixed_2),2);

for tt=1:T
    Pm_tt=mu_vary(:,tt).*conj(mu_vary(:,tt));
    search_tt=search_area2(:,tt);
    [search_tt,ind_tt]=sort(search_tt,'ascend');
    Pm_tt=Pm_tt(ind_tt);
    insert=(search_tt(1:end-1)+search_tt(2:end))/2;
    search_area_tt2=zeros(length(search_tt)*2-1,1);
    search_area_tt2(1:2:end)=search_tt;
    search_area_tt2(2:2:end)=insert;
    Pm_tt2=zeros(length(Pm_tt)*2-1,1);
    Pm_tt2(1:2:end)=Pm_tt;
    DOA(tt,:)=findmax([search_area_sort;search_area_tt2],  [Pm_fixed;Pm_tt2], N_alpha,dis);
end

%
% mu_fixed=mu(1:K_hat,:);
% Pm_fixed=sum(mu_fixed.*conj(mu_fixed),2);
% % plot(search_area,Pm_fixed)
% DOA_fixed=findmax(search_area,Pm_fixed,1,dis);
%
%
% mu_vary=mu(1+K_hat:end,:);
% for t=1:T
%     Pm_vary=sum(mu_vary(:,t).*conj(mu_vary(:,t)),2);
% %      plot(search_area2(:,t),Pm_vary)
%     DOA_vary(:,t)=findmax(search_area2(:,t),Pm_vary,1,dis);
% end
%
% DOA=[DOA_fixed*ones(T,1),DOA_vary'];









