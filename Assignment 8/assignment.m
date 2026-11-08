clc
clear all
close all

A = [-4 0; 3 -2];
C = [0 1];
Q = 400;
P0 = [0.04 0;0 0.01];
R = zeros(2,2);

%% Integrazione FDRE (punto K3)

function dPdt = riccati_eq(t, P)
    
    p11 = P(1);
    p12 = P(2);
    p22 = P(3);
    
    dp11 = -8 * p11 - 400 * p12^2;
    dp12 = 3 * p11 - 6 * p12 - 400 * p12 * p22;
    dp22 = 6 * p12 - 4 * p22 - 400 * p22^2;
    
    dPdt = [dp11; dp12; dp22];
end


tspan = [0 2];
[t, P] = ode45(@riccati_eq, tspan, [0.04; 0; 0.01]);

p11 = P(:, 1);
p12 = P(:, 2);
p22 = P(:, 3);

figure;
plot(t, p11, 'b-', t, p12, 'r--', t, p22, 'g-.', 'LineWidth', 2);
xlabel('t');
ylabel('Valori di P');
legend('p_{11}', 'p_{12}', 'p_{22}');
title('Soluzione delle equazioni di Riccati');
grid on;
