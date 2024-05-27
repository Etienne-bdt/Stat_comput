clear all
close all
clc

%% DATA SIMULATION ACCORDING TO A GMM

K = 3;
N = 1000;
d = 2;

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
hold off

%% ESTIMATION