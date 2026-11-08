clc
clear all
close all

%% Definizione sistema reale
a = 0.4;
b = 0.4;

%% Definizione Gamma
g1 = 1;
g2 = 1;

G = [g1 0; 0 g2];

%% Ritratti di fase 
% NB: questa sezione dev'essere eseguita solo dopo aver eseguito Simulink
clear theta_tilde_grad;
clear theta_tilde_DREM;

theta = [b a];

theta_tilde_grad = out.theta_hat_grad - theta;
b_tilde_grad = theta_tilde_grad(:,1);
a_tilde_grad = theta_tilde_grad(:,2);

theta_tilde_DREM = out.theta_hat_DREM - theta;
b_tilde_DREM = theta_tilde_DREM(:,1);
a_tilde_DREM = theta_tilde_DREM(:,2);

figure

x_min = -0.45;
x_max = 0.05;
y_min = -0.45;
y_max = 0.05;



subplot(2,1,1)
plot(a_tilde_grad, b_tilde_grad,'LineWidth',2)
hold on
plot(a_tilde_grad(1),b_tilde_grad(1),'o','MarkerSize',8)
plot(a_tilde_grad(end),b_tilde_grad(end),'o','MarkerSize',8)
grid on
xlabel('$\tilde a$','Interpreter','latex')
ylabel('$\tilde b$','Interpreter','latex')
title('Metodo del gradiente')
legend('Traiettoria','Inizio','Fine','Location','best')
xlim([x_min x_max])
ylim([y_min y_max])


subplot(2,1,2)
plot(a_tilde_DREM,b_tilde_DREM,'LineWidth',2)
hold on
plot(a_tilde_DREM(1),b_tilde_DREM(1),'o','MarkerSize',8)
plot(a_tilde_DREM(end),b_tilde_DREM(end),'o','MarkerSize',8)
grid on
xlabel('$\tilde a$','Interpreter','latex')
ylabel('$\tilde b$','Interpreter','latex')
title('Metodo del gradiente')
legend('Traiettoria','Inizio','Fine','Location','best')
xlim([x_min x_max])
ylim([y_min y_max])

