clc;
clear;
close all;

h= -0.1;

B=[10 -10 0 0;
    -10 20 -10 0;
    0 -10 20 -10;
    0 0 -10 10];
L = [2 -1 0 -1;
     -1 2 -1 0;
     0 -1 2 -1;
     -1 0 -1 2];
D = diag([1 2 4 8]);

A= h*B*D*L*D;



B = [ 1 1 1 1]';
C = [1 0 0 0];
D_d = [0];

sys = ss(A,B,C,D_d);

% Time vector
dt = 0.001;
t = 0:dt:10;

% Approximate impulse at t = 1
u = zeros(size(t));

% Pulse width = dt, area = 1
u(t == 1) = 1/dt;

% Simulate
[y,t,x] = lsim(sys,u,t);

x_n = x ./ [8 4 2 1];
% Plot
plot(t,x_n,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('States');
grid on;

%% Experiment 1 numerical results

% Improve and save the figure
title('A-NAPC Normalized Active-Power Response');
ylabel('Normalized Active-Power State');
legend('DER 1','DER 2','DER 3','DER 4','Location','best');
xlim([0 10]);

exportgraphics(gcf,'Experiment1_ANAPC_Baseline.png','Resolution',300);

%% 1. Calculate all system eigenvalues
lambda = eig(A);

% Remove the expected zero consensus eigenvalue
zero_tolerance = 1e-8;
lambda_effective = lambda(abs(lambda) > zero_tolerance);

% The dominant effective pole is closest to zero
[~,dominant_index] = max(real(lambda_effective));
dominant_pole = lambda_effective(dominant_index);

% Positive stability-margin convention
stability_margin = -real(dominant_pole);

%% 2. Final normalized state values
final_normalized_states = x_n(end,:);

% Difference between highest and lowest final normalized state
final_consensus_error = max(final_normalized_states) ...
                      - min(final_normalized_states);

% Mean final normalized value
final_consensus_value = mean(final_normalized_states);

%% 3. Estimate 2 percent consensus settling time
consensus_at_each_time = mean(x_n,2);

disagreement = max( ...
    abs(x_n - consensus_at_each_time), ...
    [],2);

settling_tolerance = 0.02 * abs(final_consensus_value);
start_index = find(t >= 1,1);

settling_index = NaN;

for k = start_index:length(t)
    if all(disagreement(k:end) <= settling_tolerance)
        settling_index = k;
        break
    end
end

if isnan(settling_index)
    settling_time = NaN;
else
    % Time elapsed after the disturbance at t = 1 second
    settling_time = t(settling_index) - 1;
end

%% 4. Display results in Command Window
fprintf('\nEXPERIMENT 1: A-NAPC BASELINE\n');
fprintf('--------------------------------------\n');

disp('All eigenvalues:');
disp(lambda);

disp('Effective eigenvalues:');
disp(lambda_effective);

fprintf('Dominant effective pole: %.6f %+.6fi\n', ...
    real(dominant_pole),imag(dominant_pole));

fprintf('Effective stability margin: %.6f\n', ...
    stability_margin);

fprintf('Final normalized states:\n');
disp(final_normalized_states);

fprintf('Final consensus value: %.6f\n', ...
    final_consensus_value);

fprintf('Final consensus error: %.10f\n', ...
    final_consensus_error);

fprintf('Approximate 2%% consensus settling time: %.4f s\n', ...
    settling_time);

%% 5. Create and save a results table
Metric = [
    "Stability Margin"
    "Dominant Pole Real Part"
    "Dominant Pole Imaginary Part"
    "Final Consensus Value"
    "Final Consensus Error"
    "Two Percent Settling Time"
];

Value = [
    stability_margin
    real(dominant_pole)
    imag(dominant_pole)
    final_consensus_value
    final_consensus_error
    settling_time
];

Experiment1Results = table(Metric,Value);

disp(Experiment1Results);
writetable(Experiment1Results,'Experiment1_ANAPC_Results.csv');