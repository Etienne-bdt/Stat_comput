d = 1 ; 
sigma2 = 0.5; 
N = 10000;
m=1;
m0=1;
sig0_2 = 0.5;

X = normrnd(0,sigma2,1,N);

s2 = 1/(1/(sig0_2)+N/(sigma2));
mbar = s2*((m0/sig0_2)+sum(X)/sigma2);

m_ech = normrnd(mbar,s2,1,N);

%Bizarre

%% Matrice de covariance inconnue

%p(sig|X,m) prop p(X|sig(,m))p(sig(|m))
%on trouve :
%sig suit IG((N+alpha)/2,Beta + somme des (xi-m)^2/2
%mu suit Norm((sig2m0 + Nsig02meanX)/sig2 + Nsig02 ; sig2sig02/(sig2 + N
%sig02))
alpha=2;
beta=2;


sigm = 0.5;

m_list = zeros(1,100);
sig_list = zeros(1,100);

x_g = normrnd(m, sqrt(sigma2),1,N);

for i=1:1000
       mu = normrnd((sigm * m0 + N*sig0_2*mean(x_g))/(sigm+N*sig0_2),sqrt((sigm*sig0_2)/(sigm+N*sig0_2)));
       sigm = 1./gamrnd((N+alpha)/2,2/(beta+sum((x_g-mu).^2)));
       m_list(i) = mu;
       sig_list(i) = sigm;
end

figure;
subplot(2,1,1);
plot(m_list);
title('Moyenne estimée');
subplot(2,1,2);
plot(sig_list);
title('Variance estimée');
