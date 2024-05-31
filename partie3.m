mu_0 = 1;
lambda_0 = 1;
alpha_0 = 1;
beta_0 = 1;
max_it = 50;
N = 1000;

mu = mu_0;
lambda = lambda_0;
alpha = alpha_0;
beta = beta_0;

x = normrnd(0,sqrt(0.5),N,1);


lambdal = zeros(1,max_it);
betal = zeros(1,max_it);
varl = zeros(1,max_it);
meanl = zeros(1,max_it);

for i = 1:max_it
    mu = (lambda_0*mu_0 + N*mean(x))/(lambda_0 + N);
    lambda = (lambda_0 + N)*alpha/beta;
    qmu = normrnd(mu, 1/sqrt(lambda),N,1);

    alpha = alpha_0 + N/2;
    beta = beta_0 + 0.5*sum(x.^2 - 2*x*mu + mean(qmu.^2)) + 0.5*lambda_0* (mean(qmu.^2) - 2*mu*mu_0 + mu_0^2);
    qtau = gamrnd(alpha, 1/beta,N,1);

    lambdal(i) = lambda;
    betal(i) = beta;
    varl(i) = 1/mean(qtau);
    meanl(i) = mean(qmu);
end

figure;
subplot(2,1,1);
semilogy(1:max_it,lambdal)
title('lambda');
subplot(2,1,2);
plot(betal);
title('beta');

figure;
subplot(2,1,1);
plot(meanl);
title('Moyenne estimée');
subplot(2,1,2);
plot(varl);
title('Variance estimée');
