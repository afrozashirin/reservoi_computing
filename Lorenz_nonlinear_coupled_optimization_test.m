%%This code is a variation of the code initially provided by Thomas Carroll
%%(NRL)

clear 
clc
close all

%n=100; %Nodes

y0 =  [.1 .25 1.5];

options = odeset('RelTol',1e-7,'AbsTol',1e-7);% I asked to control the ode error

[T,Y1] = ode45('Lorenz',[0:0.01: 10000],y0,options); %the time interval was [0 10000], which will create an arbitrary time step by ode45 internal time step. You need to define he time step. 
X = Y1(:,1); % Notice how I square the signal one and passed throughout the code.
Z = Y1(:,3);
meanX = mean(X); %Must change if the input signal is going to be x^2 or x^3
stdX = std(X); %Must change if the input signal is going to be x^2 or x^3
meanZ = mean(Z); %Must change if the input signal is going to be x^2 or x^3
stdZ = std(Z); %Must change if the input signal is going to be x^2 or x^3

norm_X = (X-meanX)/stdX;
norm_Y = (Y1(:,2)-mean(Y1(:,2)))/std(Y1(:,2));
norm_Z = (Z-meanZ)/stdZ;
% Tau = unifrnd(-200,200,[1 100]);
Tau = unifrnd(0,100,[1 100]);
m = max(abs(round(Tau)));

%% create the matrix for the RC
n = 100;
% linx = -100; %small
% liny = -1    ; %to large
% lambdavec  = linspace(linx,liny,n);
% A = diag(lambdavec);
Amat = load('Amat_1.mat'); 
A = Amat.A;
V = eye(n);
W=ones(n,1);
%% input signal / output signal
n_trans=60000; % transient length %40000

n_fit = 1000; % length of training signal %1000, 1200
n_test = n_fit*.5;
n_steps = n_trans+n_fit;
input_signal= norm_X;
drive_sig=input_signal(n_trans:n_steps); 
train_signal=norm_Y(n_trans:n_steps);
test_signal=norm_Y(n_steps:n_steps+n_test);

%% testing stuff
p1= 1;
p2= [0.01:0.02:5];%gamma
%p2 = 0.2;
p3= 0;
betaE = -6;

eps1= 0;% 0 for time invariant 0.5 is for time variant

for ii = 1:length(p2)
    V = eye(n);

    [T,Y] = ode45('r',[n_trans:1:n_steps]*0.01,[Y1(n_trans,:) zeros(1,n)],options,p1,p2(ii),p3,A,W,meanX,stdX,V,n); %Must change function for different attractor. 'p' is for the Hindmarsh-Rose. 'r' is for the Lorenz.
    Omega_r = zeros(n_fit, n);
    Omega_r = Y(:,4:n+3); 
    %Omega_r(:, n+1)=1; %Omega Matrix
    Omega_rdot = zeros(n_fit+1, n);
    Ydot = zeros(n_fit+1, n+3);
    for i = 1:n_fit+1
        Ydot(i,:) = r_dot(T(i),[Y(i,:)]',p1,p2(ii),p3,A,W,meanX,stdX,V,n);
    end
    Omega_rdot = Ydot(:,4:n+3); 
    
     %Omega Matrix

    M = [Omega_r  Omega_rdot ones(length(T),1)];
    
    beta = 10^(betaE);
    net_inverse=Chinv2(M,beta);
    
    coeff_vec=net_inverse*train_signal; 
    kappas = coeff_vec(1:n); %Kappa
    lambdas = coeff_vec(n+1:n+n); %Kappa
    taus = -(lambdas)./kappas;
    kappa_nplus1 = coeff_vec(2*n+1);
    %% applying timeshift
    M_delay = Omega_r-taus'.*Omega_rdot;
   % M_delay = [M_delay ones(length(T),1)];
    
    fit_signal=M_delay*kappas + kappa_nplus1*ones(length(T),1) ;  
    %fit_signal= Omega_r*kappas + Omega_rdot*lambdas + kappa_nplus1*ones(length(T),1); 

    train_err=std(fit_signal-train_signal)/std(train_signal);
    train_err_p2(ii,:) = [p2(ii) train_err]
%% Testing 
    [T_test,Y_test] = ode45('r',[n_steps:1:n_steps+n_test]*0.01,Y(end,:),options,p1,p2(ii),p3,A,W,meanX,stdX,V,n); %Must change function for different attractor. 'p' is for the Hindmarsh-Rose. 'r' is for the Lorenz.

    Omega_bar = Y_test(:,4:n+3);
    %Omega_bar(:, n+1)=1; %Omega Matrix
    
    Omega_bar_rdot = zeros(n_test+1, n);
    
    Ydot = zeros(n_test+1, n+3);
    
    for i = 1:n_test+1
        Ydot(i,:) = r_dot(T_test(i),[Y_test(i,:)]',p1,p2(ii),p3,A,W,meanX,stdX,V,n);
    end
    Omega_bar_dot = Ydot(:,4:n+3); 
    
    Mtest = [Omega_bar-taus'.*Omega_bar_dot ];
    %Pinv
    fit_signal_test =Mtest*kappas +  kappa_nplus1*ones(length(T_test),1);  
    
    %fit_signal_test= Omega_bar*kappas + Omega_bar_dot*lambdas + kappa_nplus1*ones(length(T_test),1); 
    
    test_err=std(fit_signal_test-test_signal)/std(test_signal);
    % test_errSet(i2) = test_err;
    test_err_p2(ii,:) = [p2(ii) test_err];
end

%% figure

%%
T = array2table([train_err_p2]);
Tt = array2table([test_err_p2]);

writetable(Tt,['Lorenztesting_nonlinear_coupled_optTS.dat'],'WriteVariableNames', false,'Delimiter',' ')
writetable(T,['Lorenztraining_nonlinear_coupled_optTS.dat'],'WriteVariableNames', false,'Delimiter',' ')

%writetable(Tt,['Lorenztesting_linear_optTS.dat'],'WriteVariableNames', false,'Delimiter',' ')
%writetable(T,['Lorenztraining_linear_optTS.dat'],'WriteVariableNames', false,'Delimiter',' ')


figure
hold on
plot(p2,train_err_p2(:,2),'DisplayName','Time in-variant','LineWidth',2)
set(gca, 'YScale', 'log')

xlabel('\gamma')
ylabel('\Delta_{training}')
title('Lorenz coupled Timeshift Training')
figure
plot(p2,test_err_p2(:,2),'DisplayName','Time in-variant','LineWidth',2)
set(gca, 'YScale', 'log')

xlabel('\gamma')
ylabel('\Delta_{testing}')
title('Lorenz coupled Timeshift Testing')