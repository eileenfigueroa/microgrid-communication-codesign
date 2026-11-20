%% EE 800 controlled four-bus experiments
% Extends simu_4.m without modifying Professor Wu's original file.
% Paper mapping:
%   O-NAPC, Eq. (10): A = -h * Bel * Lw * Dp
%   A-NAPC, Eq. (13): A = -h * Bel * Dp * Lw * Dp

clc;
clear;
close all;

%% Shared four-bus data
h = 0.1;  % positive controller gain; the minus sign appears in A

% Electrical-network susceptance/Laplacian matrix (chain physical grid)
Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];

% Available-capacity ratios and inverse-capacity matrix
capacity = [8 4 2 1];
Dp = diag(1 ./ capacity);

% Multiplying Dp by a common positive scalar changes time scale, not
% proportional sharing. Professor Wu's original file uses diag([1 2 4 8]),
% which equals 8*Dp. Preserve that scaling for baseline reproduction.
Dp_scaled = 8 * Dp;

% Candidate symmetric communication Laplacians
L_chain = [ 1 -1  0  0;
           -1  2 -1  0;
            0 -1  2 -1;
            0  0 -1  1];

L_cycle = [ 2 -1  0 -1;
           -1  2 -1  0;
            0 -1  2 -1;
           -1  0 -1  2];

L_complete = [ 3 -1 -1 -1;
              -1  3 -1 -1;
              -1 -1  3 -1;
              -1 -1 -1  3];

names = {'Chain','Cycle','Complete'};
Ls = {L_chain,L_cycle,L_complete};

%% Disturbance and simulation settings
dt = 0.001;
t = 0:dt:4;
u = zeros(size(t));
u(t == 1) = 1/dt;  % unit-area pulse at t = 1 s

% Same all-node disturbance channel used in the original simu_4.m
Bdist = ones(4,1);
Cout = eye(4);
Ddist = zeros(4,1);

margin_O = zeros(1,numel(Ls));
margin_A = zeros(1,numel(Ls));

%% Experiment 1: topology and control-scheme comparison
figure('Name','Normalized power sharing');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for k = 1:numel(Ls)
    Lw = Ls{k};
    A_O = -h * Bel * Lw * Dp_scaled;
    A_A = -h * Bel * Dp_scaled * Lw * Dp_scaled;

    margin_O(k) = effectiveMargin(A_O);
    margin_A(k) = effectiveMargin(A_A);

    [~,~,xO] = lsim(ss(A_O,Bdist,Cout,Ddist),u,t);
    [~,~,xA] = lsim(ss(A_A,Bdist,Cout,Ddist),u,t);
    xO_norm = xO ./ capacity;
    xA_norm = xA ./ capacity;

    nexttile(k);
    plot(t,xO_norm,'LineWidth',1.2);
    grid on;
    title(['O-NAPC  ' names{k}]);
    xlabel('Time (s)');
    ylabel('Normalized state');

    nexttile(k+3);
    plot(t,xA_norm,'LineWidth',1.2);
    grid on;
    title(['A-NAPC  ' names{k}]);
    xlabel('Time (s)');
    ylabel('Normalized state');
end
legend({'Node 1','Node 2','Node 3','Node 4'},'Location','bestoutside');

%% Experiment 2: full-system stability margin by topology
figure('Name','Full system stability margins');
topology_labels = categorical(names);
topology_labels = reordercats(topology_labels,names);

bar(topology_labels,[margin_O' margin_A']);grid on;
ylabel('Effective stability margin');
title('Four-bus closed-loop stability margin');
legend({'O-NAPC Eq. 10','A-NAPC Eq. 13'},'Location','northwest');

%% Experiment 3: gain sensitivity for the cycle topology
h_values = [0.05 0.10 0.20];
gain_margin_O = zeros(size(h_values));
gain_margin_A = zeros(size(h_values));

for k = 1:numel(h_values)
    hk = h_values(k);
    gain_margin_O(k) = effectiveMargin(-hk*Bel*L_cycle*Dp_scaled);
    gain_margin_A(k) = effectiveMargin(-hk*Bel*Dp_scaled*L_cycle*Dp_scaled);
end

figure('Name','Controller gain sensitivity');
plot(h_values,gain_margin_O,'-o','LineWidth',1.5,'MarkerSize',7);
hold on;
plot(h_values,gain_margin_A,'-s','LineWidth',1.5,'MarkerSize',7);
grid on;
xlabel('Controller gain h');
ylabel('Effective stability margin');
title('Gain sensitivity with cycle communication');
legend({'O-NAPC','A-NAPC'},'Location','northwest');

%% Numerical summary
results = table(names',margin_O',margin_A', ...
    'VariableNames',{'Topology','O_NAPC_Margin','A_NAPC_Margin'});
disp(results);

gain_results = table(h_values',gain_margin_O',gain_margin_A', ...
    'VariableNames',{'Gain_h','O_NAPC_Margin','A_NAPC_Margin'});
disp(gain_results);

%% Local function
function margin = effectiveMargin(A)
    ev = eig(A);
    effective = ev(abs(ev) > 1e-8);  % remove conservation-related zero mode
    margin = -max(real(effective));   % distance of dominant effective pole to zero
end
