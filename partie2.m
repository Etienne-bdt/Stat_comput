clear all
close all
clc

%% DATA SIMULATION ACCORDING TO A GMM

K = 3; %nombre de classe
N = 1000; %nombre d'observation 
d = 2; %taille matrice i 
id = eye(d); %matrice unité de taille d*d 

% vector of means
Vmu = [0 0 1;-1 1 0];

% vector of variances
Vsigma2 = [20 10 50]*1e-3;

% generation of the labels
Z = unidrnd(K,1,N);

% generation of the data
X = zeros(d,N);
for k=1:K
    indk = (Z==k);
    nk = sum(indk);
    X(:,indk) = Vmu(:,k) + sqrt(Vsigma2(k))*randn(d,nk);
end

% plot
Vangle = 0:0.01:2*pi; 
Tcolor = ['r', 'g', 'b'];
figure(1)
subplot(121);
hold on
for k=1:K
    indk = (Z==k);
    % data
    scatter(X(1,indk),X(2,indk),[Tcolor(k) '.'])
    % clusters
    u=Vmu(1,k)+3*sqrt(Vsigma2(k))*cos(Vangle);
    v=Vmu(2,k)+3*sqrt(Vsigma2(k))*sin(Vangle);
    plot(u,v,Tcolor(k));
end
title('Classification de référence');
hold off

%% ESTIMATION

%initialisation aléatoire des valeurs 
%initialiser avec du Kmeans (on récupere des labels, )

moy = zeros(d, K);%moy observé 
sigma = rand(1, K); %variance

% initialisation avec le Kmeans
x=zeros(2,size(v,2)) ;  
x(1,:)= u ;  
x(2,:)= v ; 
[idx,C,sumd,D] = kmeans(x',3) ; 
vect_moy = C';
vect_sig = [1,1,1] ; 
pika = [sum(idx==1)/size(u,2),sum(idx==2)/size(u,2),sum(idx==3)/size(u,2)] ; % proba qu'un pnt appartienne a une classe 

i_max = 1e5;
i=0 ;
pas_conver = 1e-16;
pas = 10 ;  %évolution a chaque itération  

vraisemblance = zeros(1, i_max) ; 
gamma = zeros(K, N); %proba a posteriori

%Algo EM
while (i<i_max && pas > pas_conver)
 i = i+1; 
 %etape E 
        for j = 1:K
        coeff = (2*pi*sigma(j))^(-d/2) * pika(j);
        expComp = exp(-0.5 / sigma(j) * sum((X - moy(:,j)).^2, 1));
        gamma(j, :) = coeff * expComp;
        end
        gamma = gamma ./ sum(gamma, 1);
% etape M, actualisation
    for j = 1:K
        Nk = sum(gamma(j, :));
        moy(:, j) = (1 / Nk) * X * gamma(j, :)'; % moy
        sigma(j) = (1 / (Nk * d)) * sum(gamma(j, :) .* sum((X - moy(:, j)).^2, 1)); % var
        pika(k) = Nk / N;
    end
      vraisemblance(i) = sum(log(sum(gamma, 1)));
      if i == 1 
           pas = 10 ; 
          else 
        pas = abs(vraisemblance(i) - vraisemblance(i - 1)); %condition de convergence  
      end
end 

disp('Moyennes théorique :');
disp(Vmu);
disp('Moyennes estimées :');
disp(moy);
disp('Variances théorique :');
disp(Vsigma2);
disp('Variances estimées :');
disp(sigma);
disp('Nombre d''itérations :')
disp(i);

subplot(122);
hold on;
Vangle = 0:0.01:2*pi;
Tcolor = ['g', 'b', 'r'];
h = zeros(3, 1); % objets graphiques légende
labels = {'Classe 1', 'Classe 2', 'Classe 3'};
for k=1:K
    indk = (Z==k);
    % data
    h(k) = scatter(X(1,indk),X(2,indk),[Tcolor(k) '.']);
    % clusters
    u = moy(1,k) + 3*sqrt(sigma(k))*cos(Vangle);
    v = moy(2,k) + 3*sqrt(sigma(k))*sin(Vangle);
    plot(u,v,Tcolor(k), 'LineWidth', 2);
end
title('Classification pour EM');
legend(h, labels, 'Location', 'best');
hold off;

% figure;
% subplot(211);
% semilogx(N_values, errors_moy, '-o');
% title('Erreur sur les moyennes');
% xlabel('Nombre d''observations (N)');
% ylabel('Erreur');
% subplot(212);
% semilogx(N_values, errors_sigma, '-o');
% title('Erreur sur les variances');
% xlabel('Nombre d''observations (N)');
% ylabel('Erreur');

%% Simu pour différents N
% Paramètres
K = 3; % classes du mélange
d = 2; % dim. des données
N_values = logspace(1, 8, 100); 

errors_moy = zeros(length(N_values), 1); % erreurs moy gaussienne
errors_sigma = zeros(length(N_values), 1); % erreurs var gaussienne

Vmu = [0 0 1;-1 1 0]; % vraies moyennes
Vsigma2 = [0.02 0.01 0.05]; % vraies variances

% Boucle différents N
for idx = 1:length(N_values)
    N = round(N_values(idx)); 
    Z = unidrnd(K,1,N); % labels de classe
    % Génération des data selon vraies moy et var
    X = zeros(d,N);
    for k=1:K
        indk = (Z==k);
        nk = sum(indk); % nbr de data dans la classe k
        X(:,indk) = Vmu(:,k) + sqrt(Vsigma2(k))*randn(d,nk); % data distrib normale
    end

    %% Estimation par algo EM

    % Initialisation avec K-means
x=zeros(2,size(v,2)) ;  
x(1,:)= u ;  
x(2,:)= v ; 
[idx,C,sumd,D] = kmeans(x',3) ; 
vect_moy = C';
vect_sig = [1,1,1] ; 
pika = [sum(idx==1)/size(u,2),sum(idx==2)/size(u,2),sum(idx==3)/size(u,2)] ; % proba qu'un pnt appartienne a une classe 

i_max = 1e5;
i=0 ;
pas_conver = 1e-16;
pas = 10 ;  %évolution a chaque itération  

vraisemblance = zeros(1, i_max) ; 
gamma = zeros(K, N); %proba a posteriori

%Algo EM
while (i<i_max && pas > pas_conver)
 i = i+1; 
 %etape E 
        for j = 1:K
        coeff = (2*pi*sigma(j))^(-d/2) * pika(j);
        expComp = exp(-0.5 / sigma(j) * sum((X - moy(:,j)).^2, 1));
        gamma(j, :) = coeff * expComp;
        end
        gamma = gamma ./ sum(gamma, 1);
% etape M, actualisation
    for j = 1:K
        Nk = sum(gamma(j, :));
        moy(:, j) = (1 / Nk) * X * gamma(j, :)'; % moy
        sigma(j) = (1 / (Nk * d)) * sum(gamma(j, :) .* sum((X - moy(:, j)).^2, 1)); % var
        pika(k) = Nk / N;
    end
      vraisemblance(i) = sum(log(sum(gamma, 1)));
      if i == 1 
           pas = 10 ; 
          else 
        pas = abs(vraisemblance(i) - vraisemblance(i - 1)); %condition de convergence  
      end
end 


    % % Affichage
    % disp('Moyennes estimées :');
    % disp(mu);
    % disp('Variances estimées :');
    % disp(sigma2);
    % disp('Poids des classes :');
    % disp(pi_k);
    % disp('Nombre d''itérations :')
    % disp(iter);

    % Calcul des erreurs
    errors_mu(idx) = norm(moy - Vmu, 'fro'); % erreur Frobenius (matrices)
    errors_sigma(idx) = norm(sigma - Vsigma2); % norme euclidienne (vecteurs)
end

% Tracé des erreurs selon N
figure;
subplot(211);
loglog(N_values, errors_mu, '-o');
title('Erreur sur les moyennes');
xlabel('Nombre d''observations (N)');
ylabel('Erreur');
subplot(212);
loglog(N_values, errors_sigma, '-o');
title('Erreur sur sigma');
ylabel('Erreur');
