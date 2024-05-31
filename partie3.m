mu_0 = 0;
lambda_0 = 0;
max_it = 1000;
N = 1000;
x = normrnd(0,sqrt(0.5),N,1);

mu = mu_0;
lambda = lambda_0;
alpha = 100;
beta = 200;


lambdal = zeros(1,max_it);
betal = zeros(1,max_it);
varl = zeros(1,max_it);

for i = 1:max_it
    mu = (lambda_0*mu_0 + N*mean(x))/(lambda_0 + N);
    lambda = (lambda_0 + N)*alpha/beta;
    qmu = normrnd(mu, 1/sqrt(lambda),N,1);

    alpha = alpha_0 + N/2;
    beta = beta_0 + 0.5*mean(sum((x - qmu).^2+lambda_0*(mu - mu_0).^2));
    qtau = gamrnd(alpha, 1/beta);

    lambdal(i) = lambda;
    betal(i) = beta;
    varl(i) = 1/mean(qtau);
    
end

figure;
subplot(2,1,1);
semilogy(1:max_it,lambdal)
title('lambda');
subplot(2,1,2);
plot(betal);
title('beta');

figure;
plot(varl);
title('variance estimée');
