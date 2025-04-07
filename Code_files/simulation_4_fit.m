clc;
clear;
close all;

% Parameters
N = 550000;           % Total population in San Francisco
latent_period = 2;    % Latent period (days)
infectious_period = 2; % Infectious period (days)
k = 1 / latent_period; % Rate of moving from E to I
g = 1 / infectious_period; % Recovery rate

% Load Data
filename = '1918_influenza.txt'; % Provide your dataset filename
data = load(filename);
time_span = data(1:20, 1);
data_cases = data(1:20, 2);

% Initial conditions
S0 = N - data_cases(1);
E0 = 0;
I0 = data_cases(1);
R0 = 0;
initial_conditions = [S0, E0, I0, R0];

% Define SEIR Model
function dydt = seir_ode(~, y, beta, k, g, N)
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    dS = -beta * S * I / N;
    dE = beta * S * I / N - k * E;
    dI = k * E - g * I;
    dR = g * I;

    dydt = [dS; dE; dI; dR];
end

% Function to solve SEIR
function cases = solve_seir(beta, time_span, initial_conditions, k, g, N)
    [~, y] = ode45(@(t, y) seir_ode(t, y, beta, k, g, N), time_span, initial_conditions);
    cases = y(:, 3); % Extract I (Infectious)
end

% Objective Function for Estimation
objective_function = @(beta, t) solve_seir(beta, t, initial_conditions, k, g, N);

% Initial guess and bounds
beta0 = 1.1;
lb = 0.1;
ub = 5.0;

% Estimating beta using lsqcurvefit
beta_estimated = lsqcurvefit(objective_function, beta0, time_span, data_cases, lb, ub);

fprintf('Estimated Transmission Rate (Beta): %.2f\n', beta_estimated);
R0_estimated = beta_estimated / g;
fprintf('Estimated R0: %.2f\n', R0_estimated);

% Bootstrap Parameters
n_bootstrap = 200;
time_days = [16, 18, 20];
R0_bootstrap = zeros(n_bootstrap, 3);
rng(100); % Reproducibility

% Perform Bootstrapping with Poisson Errors
for i = 1:3
    t_sim = time_span(1:time_days(i));
    data_sim = data_cases(1:time_days(i));

    for j = 1:n_bootstrap
        % Generate bootstrap realization using Poisson noise
        boot_data = poissrnd(data_sim);

        % Estimate beta using bootstrap data
        beta_boot = lsqcurvefit(objective_function, beta0, t_sim, boot_data, lb, ub);
        R0_bootstrap(j, i) = beta_boot / g;
    end
end

% Calculate the confidence intervals (95% CI)
lower_CI = zeros(length(time_span), 1);
upper_CI = zeros(length(time_span), 1);

boot_predictions = zeros(n_bootstrap, length(time_span));

for j = 1:n_bootstrap
    beta_boot = R0_bootstrap(j, 3) * g;
    boot_predictions(j, :) = solve_seir(beta_boot, time_span, initial_conditions, k, g, N);
end

lower_CI = prctile(boot_predictions, 2.5, 1);
upper_CI = prctile(boot_predictions, 97.5, 1);

% Plot actual vs predicted data with bootstrap confidence intervals
predicted_cases = solve_seir(beta_estimated, time_span, initial_conditions, k, g, N);

figure;
fill([time_span; flipud(time_span)], [lower_CI'; flipud(upper_CI')], 'b', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
hold on;
plot(time_span, data_cases, 'ro', 'DisplayName', 'Actual Data');
plot(time_span, predicted_cases, 'b-', 'DisplayName', 'SEIR Model');
xlabel('Days','fontweight','bold','fontsize',20);
ylabel('Infected Population','fontweight','bold','fontsize',20);
title('Actual vs Predicted Cases using SEIR Model with Bootstrap Confidence Intervals', 'FontWeight', 'bold');
legend;
grid on;
set(gca, 'FontWeight', 'bold', 'LineWidth', 2);  % Bold axis lines and ticks

% Plot the Normalized Empirical Distributions
samples_16 = R0_bootstrap(:, 1);
samples_18 = R0_bootstrap(:, 2);
samples_20 = R0_bootstrap(:, 3);
samples = {samples_16, samples_18, samples_20};
titles = {'16 Days', '18 Days', '20 Days'};

figure;
for i = 1:3
    subplot(1, 3, i);
    histogram(samples{i}, 'Normalization', 'probability', 'FaceColor', 'blue');
    xlabel('R_0','fontweight','bold','fontsize',20);
    ylabel('Probability (%)','fontweight','bold','fontsize',20);
    title(sprintf('R_0 Distribution (Time: %s)', titles{i}));
    grid on;
    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);  % Bold axis lines and ticks
end
sgtitle('Empirical Distributions of R_0 Using Different Time Periods');