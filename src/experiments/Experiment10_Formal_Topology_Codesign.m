%% EXPERIMENT 10: FORMAL A-NAPC TOPOLOGY CO-DESIGN
% EE 800 Capstone Project
%
% Extends Experiment 8 by comparing two optimization formulations:
%
% Weighted-sum form:
%   minimize_G  [1-rho(G)/rho_max] + lambda*|E(G)|/6
%
% Requirement-constrained form:
%   minimize_G  |E(G)|
%   subject to  rho(G) >= eta*rho_max,
%               T_common(G) <= T_common,max,
%               T_DER4(G) <= T_DER4,max,
%               G connected.
%
% The weighted form traces cost-performance tradeoffs. The constrained
% form answers the engineering question directly: what is the least-cost
% connected topology that satisfies explicit performance requirements?

clear; clc; close all;

%% 1. Fixed model
n = 4; h = 0.1;
capacity = [8 4 2 1];
Dp_scaled = diag([1 2 4 8]);
Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];
possible_edges = [1 2;1 3;1 4;2 3;2 4;3 4];
mmax = size(possible_edges,1);

%% 2. Disturbance and performance settings
dt = 0.001; t_end = 5; t = (0:dt:t_end)'; t_disturb = 1;
u = zeros(size(t));
[~,pulse_index] = min(abs(t-t_disturb));
u(pulse_index) = 1/dt;
B_common = ones(n,1);
B_der4 = [0;0;0;1];

%% 3. Evaluate every connected labeled four-node graph
records = struct([]); r = 0;
for mask = 1:(2^mmax-1)
    active = logical(bitget(mask,1:mmax));
    edges = possible_edges(active,:);
    Lcom = communication_laplacian(n,edges);
    if rank(Lcom) ~= n-1, continue; end

    A = -h*Bel*Dp_scaled*Lcom*Dp_scaled;
    [margin,dominant_pole] = effective_stability_margin(A);
    mc = simulate_metrics(A,B_common,capacity,t,u,t_disturb);
    md = simulate_metrics(A,B_der4,capacity,t,u,t_disturb);

    r = r+1;
    records(r).Graph_ID = r;
    records(r).Mask = mask;
    records(r).Link_Count = size(edges,1);
    records(r).Edges = edge_string(edges);
    records(r).A_NAPC_Margin = margin;
    records(r).Dominant_Real = real(dominant_pole);
    records(r).Dominant_Imag = imag(dominant_pole);
    records(r).Common_Settling_Time_s = mc.settling_time;
    records(r).DER4_Settling_Time_s = md.settling_time;
    records(r).Common_Peak_Disagreement = mc.peak_disagreement;
    records(r).DER4_Peak_Disagreement = md.peak_disagreement;
end

graphs = struct2table(records);
rho_max = max(graphs.A_NAPC_Margin);
graphs.Margin_Retained_Percent = 100*graphs.A_NAPC_Margin/rho_max;
graphs.Margin_Loss = 1-graphs.A_NAPC_Margin/rho_max;
graphs.Normalized_Link_Cost = graphs.Link_Count/mmax;

fprintf('Connected labeled graphs evaluated: %d\n',height(graphs));
fprintf('Maximum A-NAPC stability margin: %.6f\n',rho_max);

%% 4. Weighted-sum optimization (Equation 2 type)
lambda_values = [0 0.02 0.05 0.10 0.25 0.50 1.00 2.00];
weighted_records = struct([]);
for q = 1:numel(lambda_values)
    lambda = lambda_values(q);
    J = graphs.Margin_Loss+lambda*graphs.Normalized_Link_Cost;
    candidates = graphs(abs(J-min(J))<=1e-12,:);
    candidates.Objective = J(abs(J-min(J))<=1e-12);
    candidates = sortrows(candidates,{'Link_Count','A_NAPC_Margin'}, ...
        {'ascend','descend'});
    chosen = candidates(1,:);

    weighted_records(q).Lambda = lambda;
    weighted_records(q).Graph_ID = chosen.Graph_ID;
    weighted_records(q).Link_Count = chosen.Link_Count;
    weighted_records(q).Edges = chosen.Edges;
    weighted_records(q).A_NAPC_Margin = chosen.A_NAPC_Margin;
    weighted_records(q).Margin_Retained_Percent = chosen.Margin_Retained_Percent;
    weighted_records(q).Common_Settling_Time_s = chosen.Common_Settling_Time_s;
    weighted_records(q).DER4_Settling_Time_s = chosen.DER4_Settling_Time_s;
    weighted_records(q).Objective = min(J);
end
weighted = struct2table(weighted_records);

%% 5. Requirement-constrained optimization (Equation 4 type)
% Sensitivity cases show how the recommended topology changes as the
% required retained stability margin becomes stricter. The settling limits
% are explicit engineering requirements and can be revised with Dr. Wu.
eta_values = [0.90 0.93 0.94 0.95 0.98 0.99 0.999];
T_common_max = 0.40;
T_der4_max = 0.45;
constrained_records = struct([]);

for q = 1:numel(eta_values)
    eta = eta_values(q);
    feasible = graphs.A_NAPC_Margin >= eta*rho_max-1e-12 & ...
               graphs.Common_Settling_Time_s <= T_common_max+1e-12 & ...
               graphs.DER4_Settling_Time_s <= T_der4_max+1e-12;
    candidates = graphs(feasible,:);

    constrained_records(q).Eta = eta;
    constrained_records(q).Required_Margin_Percent = 100*eta;
    constrained_records(q).Common_Time_Limit_s = T_common_max;
    constrained_records(q).DER4_Time_Limit_s = T_der4_max;
    constrained_records(q).Feasible_Graph_Count = height(candidates);

    if isempty(candidates)
        constrained_records(q).Graph_ID = NaN;
        constrained_records(q).Link_Count = NaN;
        constrained_records(q).Edges = "No feasible graph";
        constrained_records(q).A_NAPC_Margin = NaN;
        constrained_records(q).Margin_Retained_Percent = NaN;
        constrained_records(q).Common_Settling_Time_s = NaN;
        constrained_records(q).DER4_Settling_Time_s = NaN;
    else
        % Primary objective: minimum links. Tie-break: maximum margin, then
        % minimum localized settling time.
        candidates = sortrows(candidates, ...
            {'Link_Count','A_NAPC_Margin','DER4_Settling_Time_s'}, ...
            {'ascend','descend','ascend'});
        chosen = candidates(1,:);
        constrained_records(q).Graph_ID = chosen.Graph_ID;
        constrained_records(q).Link_Count = chosen.Link_Count;
        constrained_records(q).Edges = chosen.Edges;
        constrained_records(q).A_NAPC_Margin = chosen.A_NAPC_Margin;
        constrained_records(q).Margin_Retained_Percent = chosen.Margin_Retained_Percent;
        constrained_records(q).Common_Settling_Time_s = chosen.Common_Settling_Time_s;
        constrained_records(q).DER4_Settling_Time_s = chosen.DER4_Settling_Time_s;
    end
end
constrained = struct2table(constrained_records);

%% 6. Recommended base case
% Base requirement: retain at least 98% of the complete-graph margin while
% meeting the two dynamic limits. This is intentionally transparent: if
% Dr. Wu prefers another requirement, change eta_base and rerun.
eta_base = 0.98;
base_row = constrained(abs(constrained.Eta-eta_base)<1e-12,:);

disp(' '); disp('WEIGHTED-SUM RESULTS');
disp(weighted(:,{'Lambda','Link_Count','Edges','Margin_Retained_Percent', ...
    'Common_Settling_Time_s','DER4_Settling_Time_s','Objective'}));
disp(' '); disp('CONSTRAINED CO-DESIGN SENSITIVITY');
disp(constrained);
disp(' '); disp('RECOMMENDED BASE CASE (eta = 0.98)');
disp(base_row);

writetable(graphs,'Experiment10_All_Graph_Performance.csv');
writetable(weighted,'Experiment10_Weighted_Selections.csv');
writetable(constrained,'Experiment10_Constrained_Selections.csv');

%% 7. Figures
% Figure 1: Constrained co-design sensitivity
fig1 = figure('Color','k','Position',[100 100 1250 620]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Left panel: selected link count
ax1 = nexttile;

stairs(ax1,100*constrained.Eta,constrained.Link_Count,'-o', ...
    'LineWidth',1.8,'MarkerSize',7);

xlabel(ax1,'Required Margin Retention, 100\eta (%)');
ylabel(ax1,'Selected Link Count');
yticks(ax1,3:6);
title(ax1,'Least-Cost Feasible Topology','Color','w');

set(ax1, ...
    'Color',[0.08 0.08 0.08], ...
    'XColor','w', ...
    'YColor','w', ...
    'GridColor',[0.55 0.55 0.55], ...
    'GridAlpha',0.35, ...
    'FontSize',11);

grid(ax1,'on');
box(ax1,'on');

if ~isempty(ax1.Toolbar)
    ax1.Toolbar.Visible = 'off';
end


% Right panel: required and achieved margin
ax2 = nexttile;

plot(ax2,100*constrained.Eta, ...
    constrained.Margin_Retained_Percent,'-s', ...
    'LineWidth',1.8,'MarkerSize',7);

hold(ax2,'on');

plot(ax2,100*constrained.Eta, ...
    100*constrained.Eta,'--', ...
    'LineWidth',1.4);

xlabel(ax2,'Required Margin Retention, 100\eta (%)');
ylabel(ax2,'Selected Margin Retention (%)');
title(ax2,'Requirement and Achieved Performance','Color','w');

legend(ax2,{'Selected topology','Minimum requirement'}, ...
    'Location','best', ...
    'TextColor','w', ...
    'Color',[0.10 0.10 0.10]);

set(ax2, ...
    'Color',[0.08 0.08 0.08], ...
    'XColor','w', ...
    'YColor','w', ...
    'GridColor',[0.55 0.55 0.55], ...
    'GridAlpha',0.35, ...
    'FontSize',11);

grid(ax2,'on');
box(ax2,'on');

if ~isempty(ax2.Toolbar)
    ax2.Toolbar.Visible = 'off';
end

title(tl,'Experiment 10: Constrained Co-Design Sensitivity', ...
    'FontWeight','bold','Color','w');

exportgraphics(fig1, ...
    'Experiment10_Constrained_Sensitivity.png', ...
    'Resolution',300, ...
    'BackgroundColor','black');

% Figure 2: Cost-performance design space
fig2 = figure('Color','k','Position',[120 120 1100 700]);
ax3 = axes(fig2);

scatter(ax3, ...
    graphs.Link_Count, ...
    graphs.Margin_Retained_Percent, ...
    75, ...
    graphs.DER4_Settling_Time_s, ...
    'filled');

hold(ax3,'on');

thresholdLine = yline(ax3,100*eta_base,'--', ...
    '98% margin requirement', ...
    'LineWidth',1.4, ...
    'Color',[0.85 0.85 0.85]);

thresholdLine.LabelHorizontalAlignment = 'right';
thresholdLine.LabelVerticalAlignment = 'bottom';

xlabel(ax3,'Communication Link Count |E(G)|');
ylabel(ax3,'A-NAPC Margin Retained (%)');

title(ax3, ...
    'Experiment 10: Feasible Cost-Performance Design Space', ...
    'Color','w');

set(ax3, ...
    'Color',[0.08 0.08 0.08], ...
    'XColor','w', ...
    'YColor','w', ...
    'GridColor',[0.55 0.55 0.55], ...
    'GridAlpha',0.35, ...
    'FontSize',11);

xlim(ax3,[2.9 6.1]);
ylim(ax3,[10 101]);
xticks(ax3,3:6);

grid(ax3,'on');
box(ax3,'on');

if ~isempty(ax3.Toolbar)
    ax3.Toolbar.Visible = 'off';
end

cb = colorbar(ax3);
cb.Color = 'w';
cb.Label.String = 'DER 4 Settling Time (s)';
cb.Label.Color = 'w';

exportgraphics(fig2, ...
    'Experiment10_Design_Space.png', ...
    'Resolution',300, ...
    'BackgroundColor','black');


disp('Files saved:');
disp('  Experiment10_All_Graph_Performance.csv');
disp('  Experiment10_Weighted_Selections.csv');
disp('  Experiment10_Constrained_Selections.csv');
disp('  Experiment10_Constrained_Sensitivity.png');
disp('  Experiment10_Design_Space.png');

%% Local functions
function L = communication_laplacian(n,edges)
    W = zeros(n);
    for q = 1:size(edges,1)
        i = edges(q,1); j = edges(q,2);
        W(i,j)=1; W(j,i)=1;
    end
    L = diag(sum(W,2))-W;
end

function text = edge_string(edges)
    labels = strings(size(edges,1),1);
    for q = 1:size(edges,1)
        labels(q)=sprintf('%d-%d',edges(q,1),edges(q,2));
    end
    text = strjoin(labels,', ');
end

function [margin,dominant_pole] = effective_stability_margin(A)
    poles = eig(A);
    [~,zero_index] = min(abs(poles));
    poles(zero_index)=[];
    [dominant_real,idx]=max(real(poles));
    dominant_pole=poles(idx); margin=-dominant_real;
end

function metrics = simulate_metrics(A,B,capacity,t,u,t_disturb)
    sys=ss(A,B,eye(size(A)),zeros(size(A,1),1));
    x=lsim(sys,u,t); xn=x./capacity;
    avg=mean(xn,2); disagreement=max(abs(xn-avg),[],2);
    final_value=mean(xn(end,:)); post=find(t>=t_disturb);
    metrics.peak_disagreement=max(disagreement(post));
    tolerance=max(0.02*abs(final_value),1e-8);
    metrics.settling_time=NaN;
    for q=post(1):numel(t)
        if all(disagreement(q:end)<=tolerance)
            metrics.settling_time=t(q)-t_disturb; break;
        end
    end
end
