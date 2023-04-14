clear;
close all;
seed_randn = 75800;    % Seed for randn and rand
randn('seed',seed_randn);
rand('seed',seed_randn);

M=16;                  % Number of sensors
T=100;                 % Number of snapshots
resolution=5;          % Grid interval
etc=10;                % Number of active grids
ex=3;                  % Number of snapshots per group
N_alpha=2;             % Number of sources
Group = floor( T/ex ); % Number of groups
T= Group*ex;    
SNR=10;                % SNR


%% Generate DOAs 
Sampling=[0:T-1];
r=0.08;   
DOA_1=-20 + 0.4*Sampling(1:T);
DOA_2= 20 - 0.4*Sampling(1:T) ;
DOA_real=[DOA_1', DOA_2'];
%% Sort DOA
index_sort=zeros(T,N_alpha);
for tt=1:T
    [DOA_real(tt,:), index_sort(tt,:) ] =   sort(DOA_real(tt,:));
end
%% Generate signal
Y=zeros(M,T); S=zeros(M,T);
for t=1:T
    Y(:,t)=signal(M, DOA_real(t,:), SNR, 1);
end
Y_norm=norm(Y,'fro')/sqrt(M*T);Y=Y/Y_norm;
Y_group=zeros(M,ex,Group);
for gg=1:Group
    Y_group(:,:,gg)= Y(:,(ex*(gg-1)+1):ex*gg) ;
end



%% Proposed
[Est_Proposed,initial]=Proposed_initial( Y_group(:,:,1),resolution,etc);
for gg=2:Group
    [temp_DOA,initial]=Proposed_tracking(Y_group(:,:,gg),resolution,etc,initial);
    Est_Proposed=[Est_Proposed;temp_DOA];
end
for tt=1:T
    [~,ind_re]=sort(index_sort(tt,:));
    Est_Proposed(tt,:)=Est_Proposed(tt,ind_re);
end
figure;
subplot(3,2,1)
plot(DOA_real,'k.');hold on;
plot(Est_Proposed,'o')
title('a) Proposed')
axis([0 100 -50 50])
ylabel('DOA')
xlabel('Snapshot')

%% 	Block AMP
Est_Block_AMP=[];
for gg=1:Group
    [temp_DOA,initial]=Proposed_initial(Y_group(:,:,gg),resolution,etc);
    Est_Block_AMP=[Est_Block_AMP;temp_DOA];
end
for tt=1:T
    [~,ind_re]=sort(index_sort(tt,:));
    Est_Block_AMP(tt,:)=Est_Block_AMP(tt,ind_re);
end
subplot(3,2,2)
plot(DOA_real,'.');hold on;
plot(Est_Block_AMP,'o')
title('b) Block AMP')
axis([0 100 -50 50])
ylabel('DOA')
xlabel('Snapshot')

%% 	Block SBL
Est_Block_SBL=[];
for gg=1:Group
    [temp_DOA]=SBL( Y_group(:,:,gg),resolution,etc,ex);
    Est_Block_SBL=[Est_Block_SBL;temp_DOA];
end
Est_Block_SBL_full=[];
for nn=1:N_alpha
    aa= repmat(Est_Block_SBL(:,nn),1,ex)';
    Est_Block_SBL_full(:,nn)=aa(:);
end
for tt=1:T
    [~,ind_re]=sort(index_sort(tt,:));
    Est_Block_SBL_full(tt,:)=Est_Block_SBL_full(tt,ind_re);
end
subplot(3,2,3)
plot(DOA_real,'.');hold on;
plot(Est_Block_SBL_full,'o')
title('c) Block SBL')
axis([0 100 -50 50])
ylabel('DOA')
xlabel('Snapshot')

%% 	Individual SBL
Est_Individual_SBL=[];
for tt=1:T
    [temp_DOA]=SBL( Y(:,tt),resolution,etc,1);
    Est_Individual_SBL=[Est_Individual_SBL;temp_DOA];
end
for tt=1:T
    [~,ind_re]=sort(index_sort(tt,:));
    Est_Individual_SBL(tt,:)=Est_Individual_SBL(tt,ind_re);
end
subplot(3,2,4)
plot(DOA_real,'.');hold on;
plot(Est_Individual_SBL,'o')
title('d) Individual SBL')
axis([0 100 -50 50])
ylabel('DOA')
xlabel('Snapshot')

%% BSPLR
[Est_BSPLR]=BSPLR(Y,resolution,etc,N_alpha);
for tt=1:T
    [~,ind_re]=sort(index_sort(tt,:));
    Est_BSPLR(tt,:)=Est_BSPLR(tt,ind_re);
end
subplot(3,2,5)
plot(DOA_real,'.');hold on;
plot(Est_BSPLR,'o')
title('e) BSPLR')
axis([0 100 -50 50])
ylabel('DOA')
xlabel('Snapshot')

%% PAST
[Est_PAST]=PAST(Y, N_alpha);
for tt=1:T
    [~,ind_re]=sort(index_sort(tt,:));
    Est_PAST(tt,:)=Est_PAST(tt,ind_re);
end
subplot(3,2,6)
plot(DOA_real,'.');hold on;
plot(Est_PAST,'o')
title('f) PAST')
axis([0 100 -50 50])
ylabel('DOA')
xlabel('Snapshot')
