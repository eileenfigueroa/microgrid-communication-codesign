%% EXPERIMENT 11: PER-NODE SETTLING-THRESHOLD VALIDATION
% EE 800 Capstone Project
%
% Purpose:
%   Follow up on Dr. Wu's feedback for Experiment 9 by replacing the
%   final-consensus-based settling threshold with a per-node disturbance
%   threshold. For a localized disturbance at DER 4, the reference is the
%   normalized nodal disturbance at DER 4. For a common disturbance, both
%   the average and maximum normalized nodal references are evaluated.
%
% Recommended stability-oriented rule:
%   Use the maximum normalized nodal disturbance for the common input and
%   the disturbed node's normalized value for the localized input.

clear; clc; close all;

%% 1. Fixed four-DER A-NAPC model
n = 4;
h = 0.1;
capacity = [8 4 2 1];
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
disturbance_vectors = {ones(n,1),[0;0;0;1]};

%% 2. Unit-area pulse setup
dt = 0.001;
t_end = 5;
t = (0:dt:t_end)';
t_disturb = 1;
u = zeros(size(t));
[~,pulse_index] = min(abs(t-t_disturb));
u(pulse_index) = 1/dt;

%% 3. Run all topologies and both disturbance scenarios
records = struct([]);
record_index = 0;

for k = 1:numel(topology_names)
    Lcom = communication_laplacian(n,edge_sets{k});
    A = -h*Bel*Dp_scaled*Lcom*Dp_scaled;
    [margin,dominant_pole] = effective_stability_margin(A);

    for s = 1:numel(scenario_names)
        B = disturbance_vectors{s};
        sys = ss(A,B,eye(n),zeros(n,1));
        x = lsim(sys,u,t);
        x_normalized = x./capacity;

        % Express the applied disturbance in the same normalized nodal
        % units used to evaluate consensus disagreement.
        normalized_nodal_disturbance = abs(B(:)./capacity(:));

        if s == 1
            % Dr. Wu's two alternatives for disturbances at all four nodes.
            average_reference = mean(normalized_nodal_disturbance);
            maximum_reference = max(normalized_nodal_disturbance);
        else
            % Only node 4 is disturbed, so use its nodal reference.
            average_reference = normalized_nodal_disturbance(4);
            maximum_reference = normalized_nodal_disturbance(4);
        end

        average_tolerance = 0.02*average_reference;
        maximum_tolerance = 0.02*maximum_reference;

        average_metrics = dynamic_metrics_with_tolerance( ...
            t,x_normalized,t_disturb,average_tolerance);
        maximum_metrics = dynamic_metrics_with_tolerance( ...
            t,x_normalized,t_disturb,maximum_tolerance);

        % Retain the original Experiment 9 metric only as a comparison.
        final_consensus_value = mean(x_normalized(end,:));
        legacy_tolerance = max(0.02*abs(final_consensus_value),1e-8);
        legacy_metrics = dynamic_metrics_with_tolerance( ...
            t,x_normalized,t_disturb,legacy_tolerance);

        record_index = record_index+1;
        records(record_index).Topology = topology_names{k};
        records(record_index).Link_Count = size(edge_sets{k},1);
        records(record_index).Scenario = scenario_names{s};
        records(record_index).Stability_Margin = margin;
        records(record_index).Dominant_Real = real(dominant_pole);
        records(record_index).Dominant_Imag = imag(dominant_pole);
        records(record_index).Final_Consensus_Value = final_consensus_value;
        records(record_index).Legacy_Tolerance = legacy_tolerance;
        records(record_index).Legacy_Settling_Time_s = ...
            legacy_metrics.settling_time;
        records(record_index).Average_Node_Reference = average_reference;
        records(record_index).Average_Node_Tolerance = average_tolerance;
        records(record_index).Average_Node_Settling_Time_s = ...
            average_metrics.settling_time;
        records(record_index).Maximum_Node_Reference = maximum_reference;
        records(record_index).Maximum_Node_Tolerance = maximum_tolerance;
        records(record_index).Maximum_Node_Settling_Time_s = ...
            maximum_metrics.settling_time;
        records(record_index).Peak_Disagreement = ...
            maximum_metrics.peak_disagreement;
    end
end

results = struct2table(records);
writetable(results,'Experiment11_Per_Node_Threshold_Results.csv');

disp(' ');
disp('EXPERIMENT 11: PER-NODE THRESHOLD RESULTS');
disp(results(:,{'Topology','Scenario','Stability_Margin', ...
    'Legacy_Tolerance','Legacy_Settling_Time_s', ...
    'Average_Node_Tolerance','Average_Node_Settling_Time_s', ...
    'Maximum_Node_Tolerance','Maximum_Node_Settling_Time_s'}));

%% 4. Recommended maximum-node comparison
common = results(strcmp(results.Scenario,scenario_names{1}),:);
local = results(strcmp(results.Scenario,scenario_names{2}),:);

recommended = table(common.Topology,common.Link_Count, ...
    common.Stability_Margin, ...
    common.Maximum_Node_Tolerance,local.Maximum_Node_Tolerance, ...
    common.Maximum_Node_Settling_Time_s, ...
    local.Maximum_Node_Settling_Time_s, ...
    local.Maximum_Node_Settling_Time_s- ...
        common.Maximum_Node_Settling_Time_s, ...
    'VariableNames',{'Topology','Link_Count','Stability_Margin', ...
    'Common_Max_Node_Tolerance','DER4_Node_Tolerance', ...
    'Common_Settling_Time_s','DER4_Settling_Time_s', ...
    'DER4_Minus_Common_s'});

writetable(recommended,'Experiment11_Recommended_Max_Node_Comparison.csv');

disp(' ');
disp('RECOMMENDED MAXIMUM-NODE THRESHOLD COMPARISON');
disp(recommended);

%% 5. Figures
fig1 = figure('Color','k','Position',[100 100 1250 650]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
bar(ax1,[common.Legacy_Settling_Time_s ...
         common.Average_Node_Settling_Time_s ...
         common.Maximum_Node_Settling_Time_s]);
set(ax1,'XTick',1:height(common),'XTickLabel',common.Topology, ...
    'XTickLabelRotation',20);
ylabel(ax1,'Consensus Settling Time (s)');
title(ax1,'Common Input: Threshold Sensitivity','Color','w');
legend(ax1,{'Legacy final-value','Average-node','Maximum-node'}, ...
    'Location','northoutside','Orientation','horizontal', ...
    'TextColor','w','Color',[0.10 0.10 0.10]);
style_dark_axes(ax1);

ax2 = nexttile;
bar(ax2,[common.Maximum_Node_Settling_Time_s ...
         local.Maximum_Node_Settling_Time_s]);
set(ax2,'XTick',1:height(common),'XTickLabel',common.Topology, ...
    'XTickLabelRotation',20);
ylabel(ax2,'Consensus Settling Time (s)');
title(ax2,'Recommended Maximum-Node Rule','Color','w');
legend(ax2,scenario_names,'Location','northoutside', ...
    'Orientation','horizontal','TextColor','w', ...
    'Color',[0.10 0.10 0.10]);
style_dark_axes(ax2);

title(tl,'Experiment 11: Per-Node Settling-Threshold Validation', ...
    'Color','w','FontWeight','bold');
exportgraphics(fig1,'Experiment11_Per_Node_Threshold_Comparison.png', ...
    'Resolution',300,'BackgroundColor','black');

fig2 = figure('Color','k','Position',[120 120 1050 650]);
ax3 = axes(fig2);
bar(ax3,[common.Maximum_Node_Tolerance local.Maximum_Node_Tolerance]);
set(ax3,'XTick',1:height(common),'XTickLabel',common.Topology, ...
    'XTickLabelRotation',20);
ylabel(ax3,'Absolute 2% Per-Node Tolerance');
title(ax3,'Experiment 11: Applied Maximum-Node Thresholds','Color','w');
legend(ax3,scenario_names,'Location','northoutside', ...
    'Orientation','horizontal','TextColor','w', ...
    'Color',[0.10 0.10 0.10]);
style_dark_axes(ax3);
exportgraphics(fig2,'Experiment11_Applied_Thresholds.png', ...
    'Resolution',300,'BackgroundColor','black');

disp('Files saved:');
disp('  Experiment11_Per_Node_Threshold_Results.csv');
disp('  Experiment11_Recommended_Max_Node_Comparison.csv');
disp('  Experiment11_Per_Node_Threshold_Comparison.png');
disp('  Experiment11_Applied_Thresholds.png');

%% Local functions
function L = communication_laplacian(n,edges)
    W = zeros(n);
    for q = 1:size(edges,1)
        i = edges(q,1); j = edges(q,2);
        W(i,j)=1; W(j,i)=1;
    end
    L = diag(sum(W,2))-W;
end

function [margin,dominant_pole] = effective_stability_margin(A)
    poles = eig(A);
    [~,zero_index] = min(abs(poles));
    poles(zero_index)=[];
    [dominant_real,index] = max(real(poles));
    dominant_pole = poles(index);
    margin = -dominant_real;
end

function metrics = dynamic_metrics_with_tolerance( ...
    t,x_normalized,t_disturb,tolerance)

    average_state = mean(x_normalized,2);
    disagreement = max(abs(x_normalized-average_state),[],2);
    post_indices = find(t>=t_disturb);

    metrics.peak_disagreement = max(disagreement(post_indices));
    metrics.settling_time = NaN;

    for q = post_indices(1):numel(t)
        if all(disagreement(q:end)<=tolerance)
            metrics.settling_time = t(q)-t_disturb;
            break;
        end
    end
end

function style_dark_axes(ax)
    ax.Color = [0.08 0.08 0.08];
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.GridColor = [0.55 0.55 0.55];
    ax.GridAlpha = 0.35;
    ax.FontSize = 11;
    grid(ax,'on');
    box(ax,'on');
    if ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end
end
