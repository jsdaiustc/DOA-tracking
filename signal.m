function [X]=signal(M, DOA_real, SNR, T)

N_alpha=length(DOA_real);
A=exp(-1i*pi*(0:M-1)'*sin(DOA_real*pi/180));
Vj=sqrt((   (10).^(SNR/10)   ) /2);
S=(randn(N_alpha,T)+1i*randn(N_alpha,T));
S=S./abs(S)*sqrt(2);
S=S*Vj;
noise=sqrt(1/2)*(randn(M,T)+1i*randn(M,T));
X=A*S+noise;
