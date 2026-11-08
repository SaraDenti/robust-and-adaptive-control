clc
clear all 
close all

s = tf('s');

beta = 2;

omega = linspace(0,20);



%% Andamento di |beta*W2*T| al variare di omega con k fissato

k_vals = [1.1 1.5 5/3 2 5];

W2 = 2 ./ (i * omega + 10);
P = 1 ./ (i * omega -1);
W1 = 1 ./ (i * omega + 1);

figure;
hold on;

for j = 1:length(k_vals)
    k = k_vals(j);
    
    L = k * P;
    T = L ./ (1+L);

    G = abs(beta .* W2 .* T);

    plot(omega,G,'LineWidth',1.8, 'DisplayName',sprintf('k = %.1f',k));
end

xlabel('\omega');
ylabel('G(\omega,k)');
legend('Location','northeast'); %la legenda deve stare in alto a destra
grid on;
hold off;


%% Grafico della funzione f(k) = 2k / (5(k-1))

% Intervallo di k 
k = linspace(1.01,5,1000);

% Funzione
f = (2*k) ./ (5*(k-1));

% Grafico
figure;
hold on;

plot(k, f, 'b', 'LineWidth', 2);
plot(k, ones(size(k)), 'r--', 'LineWidth', 2);

xlabel('k');
ylabel('f(k)');

legend('2k/(5(k-1))', 'y = 1', 'Location', 'northeast');

grid on;
xlim([1 5]);
ylim([0 5]);

hold off;

%% Diagramma di Nyquist di |beta*W2*T| per alcuni valori di k


for j = 1:length(k_vals)
    k = k_vals(j);

    % Funzione G(s)
    G = (4*k)/((s+10)*(s+k-1));

    figure

    nyquist(G)
    hold on
    grid on

    % Cerchio unitario
    theta = linspace(0,2*pi,1000);
    plot(cos(theta), sin(theta), 'r--', 'LineWidth', 1.5)

    % Valore massimo teorico del modulo
    R = max(1.2, 4*k/(10*(k-1)) + 0.2);

    axis equal
    xlim([-R R])
    ylim([-R R])
   
    title(['k = ', num2str(k)])
    legend('G(j\omega)', 'Cerchio unitario', 'Location', 'Best')
end



%% Andamento di ‖S*W1‖/(1-2‖T*W2‖ al variare di omega con k>5/3 fissato

k_vals = [2 3 4 5 6];

W2 = 2 ./ (i * omega + 10);
P = 1 ./ (i * omega -1);
W1 = 1 ./ (i * omega + 1);

figure;
hold on;

for j = 1:length(k_vals)
    k = k_vals(j);

    L = k * P;
    S = 1 ./ (1+L);
    T = L ./ (1+L);

    f1 = abs(S .* W1);
    f2 = abs(beta .* W2 .* T);

    G = f1 ./ (1-f2);

    plot(omega,G,'LineWidth',1.8, 'DisplayName',sprintf('k = %.1f',k));
end

xlabel('\omega');
ylabel('G(\omega,k)');
legend('Location','northeast'); %la legenda deve stare in alto a destra
grid on;
hold off;

%% Grafico della funzione f(k) = 5 / (3k-5)

% Intervallo di k 
k = linspace(1.7,400);

% Funzione
f = 5 ./ (3 * k - 5);

% Grafico
figure;
hold on;

plot(k, f, 'b', 'LineWidth', 2);
plot(k, ones(size(k)), 'r--', 'LineWidth', 2);

xlabel('k');
ylabel('f(k)');

legend('5/(3k-5)', 'y = 1', 'Location', 'northeast');

grid on;
xlim([1.7 400]);
ylim([0 1]);

hold off;

