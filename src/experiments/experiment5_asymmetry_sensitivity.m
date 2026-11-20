%% EE 800 Experiment 5: asymmetric communication sensitivity
% Uses the same four-bus model as Experiments 1 through 4.
% The asymmetry factor a changes the clockwise and counterclockwise
% communication weights while preserving a connected row-sum-zero cycle.

clc;
clear;
close all;

%% Fixed baseline data
h = 0.1;

% Chain electrical network
Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];

capacity = [8 4 2 1];
Dp_scaled = diag([1 2 4 8]);

% Asymmetry values. a = 0 is the symmetric cycle baseline.
a_values = [0 0.25 0.50 0.75 0.90];

%% Disturbance and simulation settings
dt = 0.001;
t = 0:dt:4;
u = zeros(size(t));
u(t == 1) = 1/dt;

Bdist = ones(4,1);
Cout = eye(4);
Ddist = zeros(4,1);

%% Storage
nA = numel(a_values);
margin_O = zeros(nA,1);
margin_A = zeros(nA,1);
dominant_real_O = zeros(nA,1);
dominant_imag_O = zeros(nA,1);
dominant_real_A = zeros(nA,1);
dominant_imag_A = zeros(nA,1);
consensus_error_O = zeros(nA,1);
consensus_error_A = zeros(nA,1);

%% Sweep the asymmetry factor
for k = 1:nA
    a = a_values(k);

    % Directed four-node cycle Laplacian.
    % Clockwise weight = 1+a; reverse weight = 1-a.
    % Every row sums to zero and a=0 recovers the symmetric cycle.
    L_asym = [2       -(1-a)   0       -(1+a);
              -(1+a)   2      -(1-a)   0;
               0       -(1+a)  2      -(1-a);
              -(1-a)   0      -(1+a)   2];

    % Check the Laplacian construction before using it.
    if max(abs(sum(L_asym,2))) > 1e-10
        error('Communication Laplacian rows do not sum to zero.');
    end

    % O-NAPC Equation 10 and A-NAPC Equation 13
    A_O = -h * Bel * L_asym * Dp_scaled;
    A_A = -h * Bel * Dp_scaled * L_asym * Dp_scaled;

    [margin_O(k),dominant_O] = effectiveMargin(A_O);
    [margin_A(k),dominant_A] = effectiveMargin(A_A);

    dominant_real_O(k) = real(dominant_O);
    dominant_imag_O(k) = imag(dominant_O);
    dominant_real_A(k) = real(dominant_A);
    dominant_imag_A(k) = imag(dominant_A);

    [~,~,xO] = lsim(ss(A_O,Bdist,Cout,Ddist),u,t);
    [~,~,xA] = lsim(ss(A_A,Bdist,Cout,Ddist),u,t);

    xO_norm = xO ./ capacity;
    xA_norm = xA ./ capacity;

    consensus_error_O(k) = max(xO_norm(end,:)) ...
                         - min(xO_norm(end,:));
    consensus_error_A(k) = max(xA_norm(end,:)) ...
                         - min(xA_norm(end,:));
end

%% Figure 1: full-system stability margin versus asymmetry
figure('Name','Asymmetry stability margins');
plot(a_values,margin_O,'-o','LineWidth',1.6,'MarkerSize',7);
hold on;
plot(a_values,margin_A,'-s','LineWidth',1.6,'MarkerSize',7);
grid on;
xlabel('Asymmetry factor a');
ylabel('Effective stability margin');
title('Stability Margin Versus Communication Asymmetry');
legend({'O-NAPC','A-NAPC'},'Location','best');
exportgraphics(gcf,'Experiment5_Asymmetry_Margins.png', ...
    'Resolution',300,'BackgroundColor','current');

%% Figure 2: representative normalized responses
selected_a = [0 0.50 0.90];
figure('Name','Asymmetric communication responses');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for k = 1:numel(selected_a)
    a = selected_a(k);
    L_asym = [2       -(1-a)   0       -(1+a);
              -(1+a)   2      -(1-a)   0;
               0       -(1+a)  2      -(1-a);
              -(1-a)   0      -(1+a)   2];

    A_O = -h * Bel * L_asym * Dp_scaled;
    A_A = -h * Bel * Dp_scaled * L_asym * Dp_scaled;

    [~,~,xO] = lsim(ss(A_O,Bdist,Cout,Ddist),u,t);
    [~,~,xA] = lsim(ss(A_A,Bdist,Cout,Ddist),u,t);

    nexttile(k);
    plot(t,xO ./ capacity,'LineWidth',1.2);
    grid on;
    xlim([0 4]);
    ylim([0 1]);
    title(sprintf('O-NAPC  a = %.2f',a));
    xlabel('Time (s)');
    ylabel('Normalized state');

    nexttile(k+3);
    plot(t,xA ./ capacity,'LineWidth',1.2);
    grid on;
    xlim([0 4]);
    ylim([0 1]);
    title(sprintf('A-NAPC  a = %.2f',a));
    xlabel('Time (s)');
    ylabel('Normalized state');
end

legend({'Node 1','Node 2','Node 3','Node 4'}, ...
    'Location','bestoutside');
exportgraphics(gcf,'Experiment5_Asymmetry_Responses.png', ...
    'Resolution',300,'BackgroundColor','current');

%% Numerical summary and CSV export
Experiment5Results = table(a_values', ...
    margin_O,margin_A, ...
    dominant_real_O,dominant_imag_O, ...
    dominant_real_A,dominant_imag_A, ...
    consensus_error_O,consensus_error_A, ...
    'VariableNames',{'Asymmetry_a', ...
    'O_NAPC_Margin','A_NAPC_Margin', ...
    'O_Dominant_Real','O_Dominant_Imag', ...
    'A_Dominant_Real','A_Dominant_Imag', ...
    'O_Final_Consensus_Error','A_Final_Consensus_Error'});

disp(Experiment5Results);
writetable(Experiment5Results,'Experiment5_Asymmetry_Results.csv');

%% Local function
function [margin,dominantPole] = effectiveMargin(A)
    ev = eig(A);
    effective = ev(abs(ev) > 1e-8);
    [~,idx] = max(real(effective));
    dominantPole = effective(idx);
    margin = -real(dominantPole);
end
