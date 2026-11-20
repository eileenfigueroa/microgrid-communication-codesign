%% EXPERIMENT 7: DYNAMIC COMPARISON OF SELECTED COMMUNICATION TOPOLOGIES
% EE 800 Capstone Project
% Compares six communication networks using the same four-DER model and
% identical disturbance scenarios. Full numerical results include O-NAPC
% and A-NAPC. Figures emphasize A-NAPC because the selected candidate
% topologies came from the Experiment 6 A-NAPC topology search.

clear; clc; close all;

%% 1. Fixed four-DER model (same model used in Experiments 1-6)
n = 4;
h = 0.1;
capacity = [8 4 2 1];
Dp_scaled = diag([1 2 4 8]);

Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];

%% 2. Selected communication topologies from Experiment 6
topology_names = { ...
    'Chain', ...
    'Symmetric Cycle', ...
    'Best 3-Link Star', ...
    'Best 4-Link Graph', ...
    'Best 5-Link Graph', ...
    'Complete Graph'};

edge_sets = { ...
    [1 2; 2 3; 3 4], ...
    [1 2; 2 3; 3 4; 1 4], ...
    [1 2; 1 3; 1 4], ...
    [1 2; 1 3; 1 4; 2 4], ...
    [1 2; 1 3; 1 4; 2 4; 3 4], ...
    [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]};

link_count = cellfun(@(E) size(E,1), edge_sets);

%% 3. Two consistent disturbance scenarios
% A unit-area pulse is applied at t = 1 s. These labels describe the input
% channels used in the mathematical state-space model; they do not assign a
% more specific physical disturbance interpretation than the model supports.
dt = 0.001;
t_end = 5;
t = (0:dt:t_end)';
t_disturb = 1;
u = zeros(size(t));
[~, pulse_index] = min(abs(t - t_disturb));
u(pulse_index) = 1/dt;

scenario_names = {'Common Input', 'Localized Input at DER 4'};
disturbance_vectors = {ones(n,1), [0; 0; 0; 1]};

controller_names = {'O-NAPC', 'A-NAPC'};

%% 4. Run every topology, controller, and disturbance combination
records = struct([]);
record_index = 0;

% Store A-NAPC trajectories for the two response figures.
a_response = cell(numel(topology_names), numel(scenario_names));

for k = 1:numel(topology_names)
    Lcom = communication_laplacian(n, edge_sets{k});

    A_O = -h * Bel * Lcom * Dp_scaled;
    A_A = -h * Bel * Dp_scaled * Lcom * Dp_scaled;
    A_models = {A_O, A_A};

    for c = 1:numel(controller_names)
        A = A_models{c};
        [margin, dominant_pole] = effective_stability_margin(A);

        for s = 1:numel(scenario_names)
            Bdist = disturbance_vectors{s};
            sys = ss(A, Bdist, eye(n), zeros(n,1));
            x = lsim(sys, u, t);
            x_normalized = x ./ capacity;

            metrics = dynamic_metrics(t, x_normalized, t_disturb);

            record_index = record_index + 1;
            records(record_index).Topology = topology_names{k};
            records(record_index).Link_Count = link_count(k);
            records(record_index).Controller = controller_names{c};
            records(record_index).Scenario = scenario_names{s};
            records(record_index).Stability_Margin = margin;
            records(record_index).Dominant_Real = real(dominant_pole);
            records(record_index).Dominant_Imag = imag(dominant_pole);
            records(record_index).Final_Consensus_Value = metrics.final_value;
            records(record_index).Final_Consensus_Error = metrics.final_error;
            records(record_index).Two_Percent_Settling_Time_s = metrics.settling_time;
            records(record_index).Peak_Disagreement = metrics.peak_disagreement;

            if strcmp(controller_names{c}, 'A-NAPC')
                a_response{k,s} = x_normalized;
            end
        end
    end
end

results = struct2table(records);
writetable(results, 'Experiment7_Dynamic_Results.csv');

%% 5. Compact A-NAPC summary table
common_settling = nan(numel(topology_names),1);
local_settling = nan(numel(topology_names),1);
common_error = nan(numel(topology_names),1);
local_error = nan(numel(topology_names),1);
common_peak = nan(numel(topology_names),1);
local_peak = nan(numel(topology_names),1);
a_margin = nan(numel(topology_names),1);

for k = 1:numel(topology_names)
    rows_common = strcmp(results.Topology, topology_names{k}) & ...
                  strcmp(results.Controller, 'A-NAPC') & ...
                  strcmp(results.Scenario, scenario_names{1});
    rows_local = strcmp(results.Topology, topology_names{k}) & ...
                 strcmp(results.Controller, 'A-NAPC') & ...
                 strcmp(results.Scenario, scenario_names{2});

    a_margin(k) = results.Stability_Margin(rows_common);
    common_settling(k) = results.Two_Percent_Settling_Time_s(rows_common);
    local_settling(k) = results.Two_Percent_Settling_Time_s(rows_local);
    common_error(k) = results.Final_Consensus_Error(rows_common);
    local_error(k) = results.Final_Consensus_Error(rows_local);
    common_peak(k) = results.Peak_Disagreement(rows_common);
    local_peak(k) = results.Peak_Disagreement(rows_local);
end

a_summary = table(topology_names', link_count', a_margin, ...
    common_settling, local_settling, common_error, local_error, ...
    common_peak, local_peak, ...
    'VariableNames', {'Topology','Link_Count','A_NAPC_Margin', ...
    'Common_Settling_Time_s','DER4_Settling_Time_s', ...
    'Common_Final_Error','DER4_Final_Error', ...
    'Common_Peak_Disagreement','DER4_Peak_Disagreement'});

writetable(a_summary, 'Experiment7_ANAPC_Summary.csv');

disp(' ');
disp('EXPERIMENT 7: A-NAPC DYNAMIC SUMMARY');
disp(a_summary);
disp('Files saved:');
disp('  Experiment7_Dynamic_Results.csv');
disp('  Experiment7_ANAPC_Summary.csv');

%% 6. Figure 1: A-NAPC response to the common input
fig1 = figure('Color','k','Position',[80 80 1450 780]);
tl1 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(tl1,'Experiment 7: A-NAPC Response to Common Unit-Area Pulse', ...
    'Color','w','FontWeight','bold');

for k = 1:numel(topology_names)
    ax = nexttile;
    plot(t, a_response{k,1}, 'LineWidth', 1.45);
    xline(t_disturb, '--', 'Pulse', 'Color',[0.9 0.9 0.9], ...
        'LabelVerticalAlignment','bottom');
    title(sprintf('%s (%d links)', topology_names{k}, link_count(k)), ...
        'Color','w');
    xlabel('Time (s)'); ylabel('Normalized State');
    ylim([-0.2 1.0]);
    legend({'DER 1','DER 2','DER 3','DER 4'}, ...
        'Location','best','TextColor','w','Color',[0.10 0.10 0.10]);
    style_dark_axes(ax);
end
exportgraphics(fig1, 'Experiment7_ANAPC_Common_Response.png', ...
    'Resolution',300,'BackgroundColor','black');

%% 7. Figure 2: A-NAPC response to localized input at DER 4
fig2 = figure('Color','k','Position',[100 100 1450 780]);
tl2 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(tl2,'Experiment 7: A-NAPC Response to Localized Unit-Area Pulse at DER 4', ...
    'Color','w','FontWeight','bold');

for k = 1:numel(topology_names)
    ax = nexttile;
    plot(t, a_response{k,2}, 'LineWidth', 1.45);
    xline(t_disturb, '--', 'Pulse', 'Color',[0.9 0.9 0.9], ...
        'LabelVerticalAlignment','bottom');
    title(sprintf('%s (%d links)', topology_names{k}, link_count(k)), ...
        'Color','w');
    xlabel('Time (s)'); ylabel('Normalized State');
    ylim([-0.2 1.0]);
    legend({'DER 1','DER 2','DER 3','DER 4'}, ...
        'Location','best','TextColor','w','Color',[0.10 0.10 0.10]);
    style_dark_axes(ax);
end
exportgraphics(fig2, 'Experiment7_ANAPC_DER4_Response.png', ...
    'Resolution',300,'BackgroundColor','black');

%% 8. Figure 3: A-NAPC settling-time comparison
fig3 = figure('Color','k','Position',[140 140 1200 650]);
ax3 = axes(fig3);
bar(ax3, [common_settling local_settling], 'grouped');
ax3.XTick = 1:numel(topology_names);
ax3.XTickLabel = topology_names;
ax3.XTickLabelRotation = 20;
ylabel(ax3,'2% Consensus Settling Time after Pulse (s)');
title(ax3,'Experiment 7: A-NAPC Dynamic Performance by Topology', ...
    'Color','w');
legend(ax3, scenario_names, 'Location','northoutside', ...
    'Orientation','horizontal','TextColor','w','Color',[0.10 0.10 0.10]);
style_dark_axes(ax3);
exportgraphics(fig3, 'Experiment7_ANAPC_Settling_Times.png', ...
    'Resolution',300,'BackgroundColor','black');

%% 9. Figure 4: communication cost versus dynamic performance
fig4 = figure('Color','k','Position',[160 160 1100 650]);
ax4 = axes(fig4);
scatter(ax4, link_count, common_settling, 105, 'o', ...
    'LineWidth',1.8); hold(ax4,'on');
scatter(ax4, link_count, local_settling, 105, 's', ...
    'LineWidth',1.8);

for k = 1:numel(topology_names)
    text(ax4, link_count(k)+0.03, local_settling(k), topology_names{k}, ...
        'Color','w','FontSize',9,'Interpreter','none');
end

xlabel(ax4,'Communication Link Count');
ylabel(ax4,'2% Consensus Settling Time after Pulse (s)');
title(ax4,'Experiment 7: Communication Cost vs. A-NAPC Dynamic Performance', ...
    'Color','w');
legend(ax4, scenario_names, 'Location','best', ...
    'TextColor','w','Color',[0.10 0.10 0.10]);
style_dark_axes(ax4);
exportgraphics(fig4, 'Experiment7_Cost_Performance_Tradeoff.png', ...
    'Resolution',300,'BackgroundColor','black');

disp('  Experiment7_ANAPC_Common_Response.png');
disp('  Experiment7_ANAPC_DER4_Response.png');
disp('  Experiment7_ANAPC_Settling_Times.png');
disp('  Experiment7_Cost_Performance_Tradeoff.png');

%% Local functions
function L = communication_laplacian(n, edges)
    W = zeros(n);
    for q = 1:size(edges,1)
        i = edges(q,1);
        j = edges(q,2);
        W(i,j) = 1;
        W(j,i) = 1;
    end
    L = diag(sum(W,2)) - W;
end

function [margin, dominant_pole] = effective_stability_margin(A)
    poles = eig(A);
    [~, zero_index] = min(abs(poles));
    effective_poles = poles;
    effective_poles(zero_index) = [];
    [dominant_real, dominant_index] = max(real(effective_poles));
    dominant_pole = effective_poles(dominant_index);
    margin = -dominant_real;
end

function metrics = dynamic_metrics(t, x_normalized, t_disturb)
    average_state = mean(x_normalized,2);
    disagreement = max(abs(x_normalized - average_state),[],2);
    final_value = mean(x_normalized(end,:));
    final_error = max(x_normalized(end,:)) - min(x_normalized(end,:));

    post_indices = find(t >= t_disturb);
    peak_disagreement = max(disagreement(post_indices));

    % The floor prevents a nearly zero final consensus value from creating
    % an unrealistically tiny numerical tolerance.
    tolerance = max(0.02 * abs(final_value), 1e-8);
    settling_time = NaN;
    for q = post_indices(1):numel(t)
        if all(disagreement(q:end) <= tolerance)
            settling_time = t(q) - t_disturb;
            break;
        end
    end

    metrics.final_value = final_value;
    metrics.final_error = final_error;
    metrics.settling_time = settling_time;
    metrics.peak_disagreement = peak_disagreement;
end

function style_dark_axes(ax)
    ax.Color = [0.08 0.08 0.08];
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.GridColor = [0.55 0.55 0.55];
    ax.MinorGridColor = [0.35 0.35 0.35];
    ax.GridAlpha = 0.35;
    ax.FontSize = 10;
    grid(ax,'on');
    box(ax,'on');
end
