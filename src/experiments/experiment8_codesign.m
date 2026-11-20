%% EXPERIMENT 8: SIMPLIFIED CONTROL-COMMUNICATION CO-DESIGN
% EE 800 Capstone Project
%
% Purpose:
%   Begin adapting the cost-performance structure of Reference [3] to the
%   validated four-DER linear A-NAPC model.
%
% Important scope statement:
%   This script does NOT reproduce the nonlinear dissipativity/LMI theorem
%   in Reference [3]. That paper co-designs controllers and communication
%   topology for a nonlinear DC microgrid with CPLs and input saturation.
%   Here, the controller model and gain remain fixed. The discrete
%   communication topology is selected by minimizing
%
%       J_lambda(G) = normalized margin loss + lambda * normalized link cost
%
%   over all 38 connected labeled four-node undirected graphs. The selected
%   designs are then validated using the two Experiment 7 disturbances.

clear; clc; close all;

%% 1. Fixed A-NAPC model
n = 4;
h = 0.1;
capacity = [8 4 2 1];
Dp_scaled = diag([1 2 4 8]);

Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];

possible_edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
number_possible_edges = size(possible_edges,1);

%% 2. Enumerate all connected labeled communication graphs
records = struct([]);
record_index = 0;

for mask = 1:(2^number_possible_edges - 1)
    active = logical(bitget(mask,1:number_possible_edges));
    edges = possible_edges(active,:);
    Lcom = communication_laplacian(n,edges);

    % A connected n-node undirected graph has Laplacian rank n-1.
    if rank(Lcom) ~= n-1
        continue;
    end

    A = -h * Bel * Dp_scaled * Lcom * Dp_scaled;
    [margin,dominant_pole] = effective_stability_margin(A);

    record_index = record_index + 1;
    records(record_index).Graph_ID = record_index;
    records(record_index).Mask = mask;
    records(record_index).Link_Count = size(edges,1);
    records(record_index).Edges = edge_string(edges);
    records(record_index).A_NAPC_Margin = margin;
    records(record_index).Dominant_Real = real(dominant_pole);
    records(record_index).Dominant_Imag = imag(dominant_pole);
end

graphs = struct2table(records);
maximum_margin = max(graphs.A_NAPC_Margin);
graphs.Margin_Loss = 1 - graphs.A_NAPC_Margin/maximum_margin;
graphs.Normalized_Link_Cost = graphs.Link_Count/number_possible_edges;

fprintf('Connected labeled graphs evaluated: %d\n',height(graphs));
fprintf('Maximum A-NAPC stability margin: %.6f\n',maximum_margin);

%% 3. Sweep the communication-cost weight lambda
% lambda = 0 selects performance only. Increasing lambda assigns greater
% importance to communication sparsity.
lambda_values = [0 0.02 0.05 0.10 0.25 0.50 1.00 2.00];
selection = struct([]);

for q = 1:numel(lambda_values)
    lambda = lambda_values(q);
    objective = graphs.Margin_Loss + ...
                lambda * graphs.Normalized_Link_Cost;

    minimum_value = min(objective);
    candidate_indices = find(abs(objective-minimum_value) <= 1e-12);

    % Deterministic tie-break: select fewer links, then larger margin.
    candidate_table = graphs(candidate_indices,:);
    candidate_table.Objective = objective(candidate_indices);
    candidate_table = sortrows(candidate_table, ...
        {'Link_Count','A_NAPC_Margin'},{'ascend','descend'});
    chosen = candidate_table(1,:);

    selection(q).Lambda = lambda;
    selection(q).Graph_ID = chosen.Graph_ID;
    selection(q).Mask = chosen.Mask;
    selection(q).Link_Count = chosen.Link_Count;
    selection(q).Edges = chosen.Edges;
    selection(q).A_NAPC_Margin = chosen.A_NAPC_Margin;
    selection(q).Margin_Retained_Percent = ...
        100*chosen.A_NAPC_Margin/maximum_margin;
    selection(q).Margin_Loss = chosen.Margin_Loss;
    selection(q).Normalized_Link_Cost = chosen.Normalized_Link_Cost;
    selection(q).Objective = minimum_value;
end

selected = struct2table(selection);

%% 4. Dynamic validation of each unique selected topology
dt = 0.001;
t_end = 5;
t = (0:dt:t_end)';
t_disturb = 1;
u = zeros(size(t));
[~,pulse_index] = min(abs(t-t_disturb));
u(pulse_index) = 1/dt;

scenario_names = {'Common Input','Localized Input at DER 4'};
disturbance_vectors = {ones(n,1),[0;0;0;1]};

selected.Common_Settling_Time_s = nan(height(selected),1);
selected.DER4_Settling_Time_s = nan(height(selected),1);
selected.Common_Peak_Disagreement = nan(height(selected),1);
selected.DER4_Peak_Disagreement = nan(height(selected),1);
selected.Common_Final_Error = nan(height(selected),1);
selected.DER4_Final_Error = nan(height(selected),1);

for q = 1:height(selected)
    active = logical(bitget(selected.Mask(q),1:number_possible_edges));
    edges = possible_edges(active,:);
    Lcom = communication_laplacian(n,edges);
    A = -h * Bel * Dp_scaled * Lcom * Dp_scaled;

    for s = 1:numel(scenario_names)
        sys = ss(A,disturbance_vectors{s},eye(n),zeros(n,1));
        x = lsim(sys,u,t);
        x_normalized = x ./ capacity;
        metrics = dynamic_metrics(t,x_normalized,t_disturb);

        if s == 1
            selected.Common_Settling_Time_s(q) = metrics.settling_time;
            selected.Common_Peak_Disagreement(q) = metrics.peak_disagreement;
            selected.Common_Final_Error(q) = metrics.final_error;
        else
            selected.DER4_Settling_Time_s(q) = metrics.settling_time;
            selected.DER4_Peak_Disagreement(q) = metrics.peak_disagreement;
            selected.DER4_Final_Error(q) = metrics.final_error;
        end
    end
end

%% 5. Save numerical results
writetable(graphs,'Experiment8_All_Graph_Codesign_Data.csv');
writetable(selected,'Experiment8_Selected_Designs.csv');

disp(' ');
disp('EXPERIMENT 8: SIMPLIFIED CO-DESIGN RESULTS');
disp(selected(:,{'Lambda','Link_Count','Edges','A_NAPC_Margin', ...
    'Margin_Retained_Percent','Objective', ...
    'Common_Settling_Time_s','DER4_Settling_Time_s'}));

%% 6. Figure 1: selected design as communication cost increases
fig1 = figure('Color','k','Position',[100 100 1250 650]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
stairs(ax1,selected.Lambda,selected.Link_Count,'-o', ...
    'LineWidth',1.8,'MarkerSize',7);
xlabel(ax1,'Communication-Cost Weight \lambda');
ylabel(ax1,'Selected Link Count');
title(ax1,'Selected Communication Cost');
yticks(ax1,3:6);
style_dark_axes(ax1);

ax2 = nexttile;
plot(ax2,selected.Lambda,selected.Margin_Retained_Percent,'-s', ...
    'LineWidth',1.8,'MarkerSize',7);
xlabel(ax2,'Communication-Cost Weight \lambda');
ylabel(ax2,'Complete-Graph Margin Retained (%)');
title(ax2,'Performance Retained by Selected Design');
ylim(ax2,[85 101]);
style_dark_axes(ax2);

title(tl,'Experiment 8: Simplified A-NAPC Communication Co-Design', ...
    'Color','w','FontWeight','bold');
exportgraphics(fig1,'Experiment8_Codesign_Selections.png', ...
    'Resolution',300,'BackgroundColor','black');

%% 7. Figure 2: objective values of the link-count champions
% For each link count, use the graph with the largest A-NAPC margin.
champions = graphs([],:);
for links = 3:6
    subset = graphs(graphs.Link_Count==links,:);
    [~,best_index] = max(subset.A_NAPC_Margin);
    champions = [champions; subset(best_index,:)]; %#ok<AGROW>
end

fig2 = figure('Color','k','Position',[120 120 1200 680]);
ax = axes(fig2);
hold(ax,'on');
line_styles = {'-o','-s','-d','-^'};

for k = 1:height(champions)
    objective_curve = champions.Margin_Loss(k) + ...
        lambda_values*champions.Normalized_Link_Cost(k);
    plot(ax,lambda_values,objective_curve,line_styles{k}, ...
        'LineWidth',1.8,'MarkerSize',7, ...
        'DisplayName',sprintf('%d-link champion',champions.Link_Count(k)));
end

xlabel(ax,'Communication-Cost Weight \lambda');
ylabel(ax,'Co-Design Objective J_\lambda');
title(ax,'Experiment 8: Performance-Cost Objective by Link Budget', ...
    'Color','w');
legend(ax,'Location','northwest','TextColor','w', ...
    'Color',[0.10 0.10 0.10]);
style_dark_axes(ax);
exportgraphics(fig2,'Experiment8_Codesign_Objective.png', ...
    'Resolution',300,'BackgroundColor','black');

disp('Files saved:');
disp('  Experiment8_All_Graph_Codesign_Data.csv');
disp('  Experiment8_Selected_Designs.csv');
disp('  Experiment8_Codesign_Selections.png');
disp('  Experiment8_Codesign_Objective.png');

%% Local functions
function L = communication_laplacian(n,edges)
    W = zeros(n);
    for q = 1:size(edges,1)
        i = edges(q,1);
        j = edges(q,2);
        W(i,j) = 1;
        W(j,i) = 1;
    end
    L = diag(sum(W,2))-W;
end

function text = edge_string(edges)
    labels = strings(size(edges,1),1);
    for q = 1:size(edges,1)
        labels(q) = sprintf('%d-%d',edges(q,1),edges(q,2));
    end
    text = strjoin(labels,', ');
end

function [margin,dominant_pole] = effective_stability_margin(A)
    poles = eig(A);
    [~,zero_index] = min(abs(poles));
    effective_poles = poles;
    effective_poles(zero_index) = [];
    [dominant_real,dominant_index] = max(real(effective_poles));
    dominant_pole = effective_poles(dominant_index);
    margin = -dominant_real;
end

function metrics = dynamic_metrics(t,x_normalized,t_disturb)
    average_state = mean(x_normalized,2);
    disagreement = max(abs(x_normalized-average_state),[],2);
    final_value = mean(x_normalized(end,:));
    final_error = max(x_normalized(end,:))-min(x_normalized(end,:));
    post_indices = find(t>=t_disturb);
    peak_disagreement = max(disagreement(post_indices));
    tolerance = max(0.02*abs(final_value),1e-8);
    settling_time = NaN;
    for q = post_indices(1):numel(t)
        if all(disagreement(q:end)<=tolerance)
            settling_time = t(q)-t_disturb;
            break;
        end
    end
    metrics.final_error = final_error;
    metrics.settling_time = settling_time;
    metrics.peak_disagreement = peak_disagreement;
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
end
