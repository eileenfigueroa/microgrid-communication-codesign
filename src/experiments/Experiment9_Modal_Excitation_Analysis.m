%% EXPERIMENT 9: MODAL EXCITATION ANALYSIS OF COMMON AND LOCALIZED INPUTS
% EE 800 Capstone Project
%
% Purpose:
%   Explain why the common input in Experiment 7 can produce a shorter
%   measured consensus settling time than the localized DER 4 input even
%   though both tests use the same state matrix and therefore have the same
%   eigenvalues and eigenvalue-based stability margin.
%
% Main idea:
%   Eigenvalues describe the decay rates available in the system. The input
%   vector determines how strongly each mode is excited. The output metric
%   (normalized-state disagreement) determines how strongly each mode is
%   observed. Thus, identical poles do not imply identical transients.

clear; clc; close all;

%% 1. Fixed model and selected A-NAPC topologies from Experiment 7
n = 4;
h = 0.1;
capacity = [8 4 2 1];
Cnorm = diag(1./capacity);
Dp_scaled = diag([1 2 4 8]);

Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];

topology_names = {'Chain','Symmetric Cycle','Best 3-Link Star', ...
    'Best 4-Link Graph','Best 5-Link Graph','Complete Graph'};

edge_sets = { ...
    [1 2; 2 3; 3 4], ...
    [1 2; 2 3; 3 4; 1 4], ...
    [1 2; 1 3; 1 4], ...
    [1 2; 1 3; 1 4; 2 4], ...
    [1 2; 1 3; 1 4; 2 4; 3 4], ...
    [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]};

scenario_names = {'Common Input','Localized Input at DER 4'};
Bset = {ones(n,1),[0;0;0;1]};

% P removes the average normalized state, so z = P*Cnorm*x is the
% disagreement output used to assess consensus.
P = eye(n) - ones(n)/n;
Cdis = P*Cnorm;

%% 2. Time-domain setup (same unit-area pulse as Experiment 7)
dt = 0.001;
t_end = 5;
t = (0:dt:t_end)';
t_disturb = 1;
u = zeros(size(t));
[~,pulse_index] = min(abs(t-t_disturb));
u(pulse_index) = 1/dt;

records = struct([]);
record_index = 0;
response = cell(numel(topology_names),numel(scenario_names));

%% 3. Eigenvalue, modal-excitation, and dynamic analysis
for k = 1:numel(topology_names)
    Lcom = communication_laplacian(n,edge_sets{k});
    A = -h*Bel*Dp_scaled*Lcom*Dp_scaled;

    % Right eigenvectors V and modal coordinates alpha = inv(V)*B.
    [V,D] = eig(A);
    poles = diag(D);
    Vinv = inv(V);

    [~,zero_index] = min(abs(poles));
    effective_indices = setdiff(1:n,zero_index);
    [~,order] = sort(real(poles(effective_indices)),'descend');
    effective_indices = effective_indices(order);
    dominant_index = effective_indices(1);
    margin = -real(poles(dominant_index));

    for s = 1:numel(scenario_names)
        B = Bset{s};
        alpha = Vinv*B;

        % A mode must be both excited by B and visible in disagreement.
        mode_visibility = vecnorm(Cdis*V,2,1).';
        contribution = abs(alpha).*mode_visibility;
        effective_contribution = contribution(effective_indices);
        total_effective = sum(effective_contribution);
        if total_effective > 0
            effective_fraction = effective_contribution/total_effective;
        else
            effective_fraction = zeros(size(effective_contribution));
        end

        dominant_excitation = abs(alpha(dominant_index));
        dominant_visibility = mode_visibility(dominant_index);
        dominant_contribution = contribution(dominant_index);
        dominant_fraction = effective_fraction(1);

        sys = ss(A,B,eye(n),zeros(n,1));
        x = lsim(sys,u,t);
        x_normalized = x./capacity;
        metrics = dynamic_metrics(t,x_normalized,t_disturb);
        response{k,s} = x_normalized;

        record_index = record_index+1;
        records(record_index).Topology = topology_names{k};
        records(record_index).Link_Count = size(edge_sets{k},1);
        records(record_index).Scenario = scenario_names{s};
        records(record_index).Stability_Margin = margin;
        records(record_index).Dominant_Pole_Real = real(poles(dominant_index));
        records(record_index).Dominant_Pole_Imag = imag(poles(dominant_index));
        records(record_index).Dominant_Mode_Input_Projection = dominant_excitation;
        records(record_index).Dominant_Mode_Output_Visibility = dominant_visibility;
        records(record_index).Dominant_Mode_Contribution = dominant_contribution;
        records(record_index).Dominant_Mode_Fraction_Percent = 100*dominant_fraction;
        records(record_index).Final_Consensus_Value = metrics.final_value;
        records(record_index).Absolute_Two_Percent_Tolerance = metrics.tolerance;
        records(record_index).Settling_Time_s = metrics.settling_time;
        records(record_index).Peak_Normalized_Settling_Time_s = ...
            metrics.peak_normalized_settling_time;
        records(record_index).Peak_Disagreement = metrics.peak_disagreement;
    end
end

results = struct2table(records);
writetable(results,'Experiment9_Modal_Summary.csv');

disp(' ');
disp('EXPERIMENT 9: MODAL EXCITATION SUMMARY');
disp(results(:,{'Topology','Scenario','Stability_Margin', ...
    'Dominant_Mode_Input_Projection','Dominant_Mode_Output_Visibility', ...
    'Dominant_Mode_Contribution','Dominant_Mode_Fraction_Percent', ...
    'Final_Consensus_Value','Absolute_Two_Percent_Tolerance', ...
    'Settling_Time_s','Peak_Normalized_Settling_Time_s','Peak_Disagreement'}));

%% 4. Paired comparison used to answer Dr. Wu's question
common_rows = strcmp(results.Scenario,scenario_names{1});
local_rows = strcmp(results.Scenario,scenario_names{2});
common = results(common_rows,:);
local = results(local_rows,:);

comparison = table(common.Topology,common.Link_Count,common.Stability_Margin, ...
    common.Dominant_Mode_Contribution,local.Dominant_Mode_Contribution, ...
    local.Dominant_Mode_Contribution./max(common.Dominant_Mode_Contribution,eps), ...
    common.Absolute_Two_Percent_Tolerance,local.Absolute_Two_Percent_Tolerance, ...
    common.Settling_Time_s,local.Settling_Time_s, ...
    common.Peak_Normalized_Settling_Time_s, ...
    local.Peak_Normalized_Settling_Time_s, ...
    local.Settling_Time_s-common.Settling_Time_s, ...
    'VariableNames',{'Topology','Link_Count','Stability_Margin', ...
    'Common_Dominant_Contribution','DER4_Dominant_Contribution', ...
    'DER4_to_Common_Contribution_Ratio','Common_Absolute_Tolerance', ...
    'DER4_Absolute_Tolerance','Common_Settling_Time_s', ...
    'DER4_Settling_Time_s','Common_Peak_Normalized_Settling_Time_s', ...
    'DER4_Peak_Normalized_Settling_Time_s', ...
    'DER4_Settling_Time_Increase_s'});
writetable(comparison,'Experiment9_Common_vs_DER4.csv');

disp(' ');
disp('PAIRED COMMON-VERSUS-DER4 COMPARISON');
disp(comparison);

%% 5. Figures
% Figure 1: Modal contribution and settling time
fig1 = figure('Color','k','Position',[100 100 1250 620]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
bar(ax1,[common.Dominant_Mode_Contribution ...
         local.Dominant_Mode_Contribution]);

set(ax1,'XTick',1:height(common), ...
    'XTickLabel',common.Topology, ...
    'XTickLabelRotation',20, ...
    'Color',[0.08 0.08 0.08], ...
    'XColor','w','YColor','w');

ylabel(ax1,'Dominant Disagreement-Mode Contribution');
title(ax1,'Input Projection and Output Visibility','Color','w');
legend(ax1,scenario_names,'Location','northoutside', ...
    'Orientation','horizontal','TextColor','w', ...
    'Color',[0.10 0.10 0.10]);
grid(ax1,'on');
box(ax1,'on');

ax2 = nexttile;
bar(ax2,[common.Settling_Time_s local.Settling_Time_s]);

set(ax2,'XTick',1:height(common), ...
    'XTickLabel',common.Topology, ...
    'XTickLabelRotation',20, ...
    'Color',[0.08 0.08 0.08], ...
    'XColor','w','YColor','w');

ylabel(ax2,'2% Consensus Settling Time (s)');
title(ax2,'Measured Dynamic Response','Color','w');
legend(ax2,scenario_names,'Location','northoutside', ...
    'Orientation','horizontal','TextColor','w', ...
    'Color',[0.10 0.10 0.10]);
grid(ax2,'on');
box(ax2,'on');

title(tl,'Experiment 9: Same Poles, Different Modal Excitation', ...
    'FontWeight','bold','Color','w');

exportgraphics(fig1, ...
    'Experiment9_Modal_Excitation_Comparison.png', ...
    'Resolution',300,'BackgroundColor','black');

% Figure 2: Modal contribution versus settling time
fig2 = figure('Color','k','Position',[120 120 1000 650]);
ax3 = axes(fig2);

scatter(ax3,common.Dominant_Mode_Contribution, ...
    common.Settling_Time_s,95,'o','LineWidth',1.7);
hold(ax3,'on');

scatter(ax3,local.Dominant_Mode_Contribution, ...
    local.Settling_Time_s,95,'s','LineWidth',1.7);

% Use numbers instead of long, overlapping topology names
for k = 1:height(common)
    text(ax3,common.Dominant_Mode_Contribution(k), ...
        common.Settling_Time_s(k),sprintf('  %d',k), ...
        'Color','w','FontSize',10,'FontWeight','bold');

    text(ax3,local.Dominant_Mode_Contribution(k), ...
        local.Settling_Time_s(k),sprintf('  %d',k), ...
        'Color','w','FontSize',10,'FontWeight','bold');
end

xlabel(ax3,'Dominant Disagreement-Mode Contribution');
ylabel(ax3,'2% Consensus Settling Time (s)');
title(ax3,'Experiment 9: Modal Excitation versus Settling Time', ...
    'Color','w');

legend(ax3,scenario_names,'Location','best', ...
    'TextColor','w','Color',[0.10 0.10 0.10]);

set(ax3,'Color',[0.08 0.08 0.08], ...
    'XColor','w','YColor','w');

grid(ax3,'on');
box(ax3,'on');

% Topology number key
topology_key = {
    '1 = Chain'
    '2 = Symmetric Cycle'
    '3 = Best 3-Link Star'
    '4 = Best 4-Link Graph'
    '5 = Best 5-Link Graph'
    '6 = Complete Graph'};

annotation(fig2,'textbox',[0.66 0.67 0.25 0.19], ...
    'String',topology_key, ...
    'FitBoxToText','on', ...
    'Color','w', ...
    'BackgroundColor',[0.10 0.10 0.10], ...
    'EdgeColor',[0.70 0.70 0.70], ...
    'FontSize',9);

exportgraphics(fig2, ...
    'Experiment9_Contribution_vs_Settling.png', ...
    'Resolution',300,'BackgroundColor','black');


disp('Files saved:');
disp('  Experiment9_Modal_Summary.csv');
disp('  Experiment9_Common_vs_DER4.csv');
disp('  Experiment9_Modal_Excitation_Comparison.png');
disp('  Experiment9_Contribution_vs_Settling.png');

%% Local functions
function L = communication_laplacian(n,edges)
    W = zeros(n);
    for q = 1:size(edges,1)
        i = edges(q,1); j = edges(q,2);
        W(i,j) = 1; W(j,i) = 1;
    end
    L = diag(sum(W,2))-W;
end

function metrics = dynamic_metrics(t,x_normalized,t_disturb)
    average_state = mean(x_normalized,2);
    disagreement = max(abs(x_normalized-average_state),[],2);
    final_value = mean(x_normalized(end,:));
    post_indices = find(t>=t_disturb);
    metrics.peak_disagreement = max(disagreement(post_indices));
    tolerance = max(0.02*abs(final_value),1e-8);
    peak_tolerance = max(0.02*metrics.peak_disagreement,1e-8);
    metrics.final_value = final_value;
    metrics.tolerance = tolerance;
    metrics.settling_time = NaN;
    metrics.peak_normalized_settling_time = NaN;
    for q = post_indices(1):numel(t)
        if all(disagreement(q:end)<=tolerance)
            metrics.settling_time = t(q)-t_disturb;
            break;
        end
    end
    for q = post_indices(1):numel(t)
        if all(disagreement(q:end)<=peak_tolerance)
            metrics.peak_normalized_settling_time = t(q)-t_disturb;
            break;
        end
    end
end
