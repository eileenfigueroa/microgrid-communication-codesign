%% EXPERIMENT 12: JOINT COMMUNICATION-LINK AND WEIGHT OPTIMIZATION
% EE 800 Capstone Project
%
% Dr. Wu's clarification is implemented directly:
%   1) The A-NAPC controller gain and physical-system setup are fixed.
%   2) The six undirected communication weights are decision variables.
%   3) A zero (or numerically negligible) weight means the link is absent.
%   4) The selected design must satisfy stability, settling-time, and
%      connectivity requirements.
%
% Decision vector:
%   w = [w12 w13 w14 w23 w24 w34]',  0 <= wij <= 1.
%
% Sparse constrained formulation (smooth zero-norm approximation):
%   minimize_w  sum_e (1-exp(-w_e/tau)) + epsilon*sum_e w_e
%   subject to  rho(w) >= eta*rho_ref,
%               T_common(w) <= T_common,max,
%               T_DER4(w) <= T_DER4,max,
%               lambda_2(L(w)) >= connectivity_min.
%
% The first objective term approximates the number of nonzero links. The
% small second term breaks ties in favor of lower total communication
% weight. Per-node 2% tolerances follow Experiment 11.

clear; clc; close all;
rng(12,'twister');

%% 1. Fixed A-NAPC model and optimization settings
model.n = 4;
model.h = 0.1;
model.capacity = [8 4 2 1];
model.Dp_scaled = diag([1 2 4 8]);
model.Bel = [ 10 -10   0   0;
             -10  20 -10   0;
               0 -10  20 -10;
               0   0 -10  10];
model.possible_edges = [1 2;1 3;1 4;2 3;2 4;3 4];
model.m = size(model.possible_edges,1);

% Controller and physical-system parameters above remain fixed.
% Only communication weights are optimized.
lb = zeros(model.m,1);
ub = ones(model.m,1);

% Engineering requirements used in Experiment 10.
eta = 0.98;
T_common_max = 0.40;
T_der4_max = 0.45;
connectivity_min = 1e-3;

% Experiment 11 maximum-node threshold rule.
model.B_common = ones(model.n,1);
model.B_der4 = [0;0;0;1];
model.common_tolerance = 0.02*max(abs(model.B_common(:)./model.capacity(:)));
model.der4_tolerance = 0.02*abs(model.B_der4(4)/model.capacity(4));

% Fast optimization grid; the final design is validated on a finer grid.
model.t_opt = (0:0.005:3)';
model.t_validate = (0:0.001:5)';

% Complete unit-weight graph defines the reference stability margin.
w_complete = ones(model.m,1);
ref = evaluate_design(w_complete,model,model.t_validate);
rho_ref = ref.margin;
rho_required = eta*rho_ref;

requirements = struct('rho_required',rho_required, ...
    'T_common_max',T_common_max,'T_der4_max',T_der4_max, ...
    'connectivity_min',connectivity_min);

fprintf('Reference complete-graph margin: %.6f\n',rho_ref);
fprintf('Required margin (eta = %.3f): %.6f\n',eta,rho_required);
fprintf('Per-node tolerances: common = %.4f, DER 4 = %.4f\n', ...
    model.common_tolerance,model.der4_tolerance);

%% 2. Sparse optimization over all six link weights
% Continuation from a smoother to a sharper approximation of link count.
tau_values = [0.20 0.10 0.05 0.02];
weight_tie_break = 0.01;

% Deterministic starts include the known Experiment 10 designs plus random
% weighted graphs. These are starting points, not enumerated solutions.
starts = [ ...
    1 1 1 1 1 1; ...                         % complete graph
    1 1 1 0 1 1; ...                         % best 5-link graph
    1 1 1 0 1 0; ...                         % best 4-link graph
    1 1 1 0 0 0; ...                         % 3-link star
    1 0 0 1 0 1; ...                         % chain
    0.65+0.35*rand(3,model.m)]';

use_patternsearch = exist('patternsearch','file') == 2;
if use_patternsearch
    fprintf('Solver: patternsearch (recommended derivative-free solver)\n');
    solver_options = optimoptions('patternsearch', ...
        'Display','off','MaxIterations',180, ...
        'MaxFunctionEvaluations',3500,'MeshTolerance',1e-5);
else
    fprintf(['Solver: fmincon fallback. Install Global Optimization ', ...
        'Toolbox for patternsearch if convergence is poor.\n']);
    solver_options = optimoptions('fmincon','Algorithm','sqp', ...
        'Display','off','MaxIterations',350, ...
        'MaxFunctionEvaluations',8000, ...
        'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-6, ...
        'StepTolerance',1e-9);
end

candidate_records = struct([]);
candidate_index = 0;

for s = 1:size(starts,2)
    w_current = min(max(starts(:,s),lb),ub);

    for q = 1:numel(tau_values)
        tau = tau_values(q);
        objective = @(w) sparse_objective(w,tau,weight_tie_break);
        nonlcon = @(w) performance_constraints(w,model,requirements,model.t_opt);

        if use_patternsearch
            [w_current,~,exitflag] = patternsearch(objective,w_current, ...
                [],[],[],[],lb,ub,nonlcon,solver_options);
        else
            [w_current,~,exitflag] = fmincon(objective,w_current, ...
                [],[],[],[],lb,ub,nonlcon,solver_options);
        end
    end

    metrics = evaluate_design(w_current,model,model.t_validate);
    [c,~] = performance_constraints(w_current,model,requirements,model.t_validate);

    candidate_index = candidate_index+1;
    candidate_records(candidate_index).Start_ID = s;
    candidate_records(candidate_index).Exit_Flag = exitflag;
    candidate_records(candidate_index).Feasible = all(c <= 1e-6);
    candidate_records(candidate_index).Smooth_Link_Cost = ...
        sum(1-exp(-w_current/tau_values(end)));
    candidate_records(candidate_index).Total_Weight = sum(w_current);
    candidate_records(candidate_index).Margin = metrics.margin;
    candidate_records(candidate_index).Margin_Retained_Percent = ...
        100*metrics.margin/rho_ref;
    candidate_records(candidate_index).Common_Settling_Time_s = ...
        metrics.common_settling_time;
    candidate_records(candidate_index).DER4_Settling_Time_s = ...
        metrics.der4_settling_time;
    candidate_records(candidate_index).Algebraic_Connectivity = ...
        metrics.algebraic_connectivity;
    candidate_records(candidate_index).Weights = w_current(:)';
end

candidates = struct2table(candidate_records);
feasible_rows = find(candidates.Feasible);
if isempty(feasible_rows)
    error(['No feasible weighted design was found. Try increasing ', ...
        'MaxFunctionEvaluations or relaxing the requirements.']);
end

% Select the smallest smooth link cost; use total weight and margin as
% transparent tie-breakers.
feasible_table = candidates(feasible_rows,:);
feasible_table = sortrows(feasible_table, ...
    {'Smooth_Link_Cost','Total_Weight','Margin'}, ...
    {'ascend','ascend','descend'});
best_raw_weights = feasible_table.Weights(1,:)';

%% 3. Convert negligible weights to absent links and refine
% The threshold is numerical, not a physical minimum communication weight.
active_threshold = 0.02;
active = best_raw_weights > active_threshold;
thresholded_weights = best_raw_weights;
thresholded_weights(~active) = 0;

% If hard thresholding breaks feasibility, retain links in descending
% weight order until all requirements are restored.
[c_threshold,~] = performance_constraints(thresholded_weights,model, ...
    requirements,model.t_validate);
if any(c_threshold > 1e-6)
    [~,order] = sort(best_raw_weights,'descend');
    thresholded_weights(:) = 0;
    for k = 1:model.m
        thresholded_weights(order(k)) = best_raw_weights(order(k));
        [c_threshold,~] = performance_constraints(thresholded_weights, ...
            model,requirements,model.t_validate);
        if all(c_threshold <= 1e-6), break; end
    end
end

active = thresholded_weights > 0;
refine_lb = zeros(model.m,1);
refine_ub = ones(model.m,1);
refine_ub(~active) = 0;
refine_objective = @(w) sum(w);
refine_nonlcon = @(w) performance_constraints(w,model,requirements,model.t_opt);

if use_patternsearch
    [w_optimized,~,refine_exitflag] = patternsearch(refine_objective, ...
        thresholded_weights,[],[],[],[],refine_lb,refine_ub, ...
        refine_nonlcon,solver_options);
else
    [w_optimized,~,refine_exitflag] = fmincon(refine_objective, ...
        thresholded_weights,[],[],[],[],refine_lb,refine_ub, ...
        refine_nonlcon,solver_options);
end

% Guard against small validation-grid violations after coarse optimization.
optimized = evaluate_design(w_optimized,model,model.t_validate);
[c_final,~] = performance_constraints(w_optimized,model,requirements, ...
    model.t_validate);
if any(c_final > 1e-5)
    warning(['Refined result did not pass the fine-grid validation; ', ...
        'reporting the feasible pre-refinement design instead.']);
    w_optimized = thresholded_weights;
    optimized = evaluate_design(w_optimized,model,model.t_validate);
    refine_exitflag = NaN;
end

% Values below numerical_zero are solver noise and are reported as zero.
% Revalidate after this final numerical cleanup.
numerical_zero = 1e-6;
active = w_optimized > numerical_zero;
w_optimized(~active) = 0;
optimized = evaluate_design(w_optimized,model,model.t_validate);
[c_cleanup,~] = performance_constraints(w_optimized,model,requirements, ...
    model.t_validate);
if any(c_cleanup > 1e-5)
    warning(['Numerical zero cleanup affected feasibility; restoring the ', ...
        'feasible thresholded design.']);
    w_optimized = thresholded_weights;
    active = w_optimized > numerical_zero;
    optimized = evaluate_design(w_optimized,model,model.t_validate);
end

%% 4. Compare with Experiment 10 fixed unit-weight base design
w_exp10 = [1;1;1;0;1;0]; % edges 1-2, 1-3, 1-4, 2-4
exp10 = evaluate_design(w_exp10,model,model.t_validate);

design_names = ["Experiment 10 unit weights";"Experiment 12 optimized weights"];
comparison = table(design_names, ...
    [nnz(w_exp10);nnz(active)], ...
    [sum(w_exp10);sum(w_optimized)], ...
    [exp10.margin;optimized.margin], ...
    100*[exp10.margin;optimized.margin]/rho_ref, ...
    [exp10.common_settling_time;optimized.common_settling_time], ...
    [exp10.der4_settling_time;optimized.der4_settling_time], ...
    [exp10.algebraic_connectivity;optimized.algebraic_connectivity], ...
    'VariableNames',{'Design','Active_Link_Count','Total_Weight', ...
    'Stability_Margin','Margin_Retained_Percent', ...
    'Common_Settling_Time_s','DER4_Settling_Time_s', ...
    'Algebraic_Connectivity'});

edge_labels = strings(model.m,1);
for k = 1:model.m
    edge_labels(k) = sprintf('%d-%d',model.possible_edges(k,1), ...
        model.possible_edges(k,2));
end

weight_table = table(edge_labels,w_exp10,w_optimized,active, ...
    'VariableNames',{'Edge','Experiment10_Weight','Optimized_Weight', ...
    'Selected_Link'});

disp(' ');
disp('EXPERIMENT 12: OPTIMIZED COMMUNICATION DESIGN');
disp(weight_table);
disp(' ');
disp('COMPARISON WITH EXPERIMENT 10');
disp(comparison);
fprintf('Refinement exit flag: %g\n',refine_exitflag);

writetable(candidates,'Experiment12_Optimization_Starts.csv');
writetable(weight_table,'Experiment12_Optimized_Weights.csv');
writetable(comparison,'Experiment12_Design_Comparison.csv');

%% 5. Figures
% Publication-ready colors. The solver and numerical results above are not
% changed; only the presentation of the two exported figures is updated.
blue   = [0.000 0.443 0.737];
orange = [0.850 0.325 0.098];
maroon = [0.650 0.000 0.200];
gold   = [1.000 0.760 0.100];
gray   = [0.350 0.350 0.350];

% Figure 1: optimized topology with edge widths proportional to weights.
fig1 = figure('Color','w','Position',[100 100 900 720]);
ax1 = axes(fig1,'Position',[0.08 0.08 0.84 0.78]);
hold(ax1,'on'); axis(ax1,'equal'); axis(ax1,'off');
ax1.Color = 'w';
if ~isempty(ax1.Toolbar), ax1.Toolbar.Visible = 'off'; end

% Diamond arrangement prevents the selected links from crossing:
% DER 1 at top, DER 2 at left, DER 3 at right, and DER 4 at bottom.
node_xy = [0.5 1;0 0.5;1 0.5;0.5 0];
for e = 1:model.m
    if w_optimized(e) <= numerical_zero, continue; end
    i = model.possible_edges(e,1); j = model.possible_edges(e,2);
    plot(ax1,node_xy([i j],1),node_xy([i j],2),'-', ...
        'Color',maroon,'LineWidth',2.5+5.5*w_optimized(e));
    edge_vector = node_xy(j,:)-node_xy(i,:);
    edge_normal = [-edge_vector(2) edge_vector(1)]/norm(edge_vector);
    label_position = mean(node_xy([i j],:),1)+0.050*edge_normal;
    text(ax1,label_position(1),label_position(2), ...
        sprintf('%.3f',w_optimized(e)), ...
        'Color','k','FontSize',12,'FontWeight','bold', ...
        'BackgroundColor','w','Margin',3, ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end
scatter(ax1,node_xy(:,1),node_xy(:,2),1500,gold, ...
    'filled','MarkerEdgeColor',maroon,'LineWidth',2.2);
for i = 1:model.n
    text(ax1,node_xy(i,1),node_xy(i,2),sprintf('DER %d',i), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontWeight','bold','FontSize',10.5,'Color','k');
end
xlim(ax1,[-0.18 1.18]); ylim(ax1,[-0.14 1.14]);
title(ax1,'Experiment 12: Optimized Communication Topology', ...
    'Color','k','FontWeight','bold','FontSize',17);
subtitle(ax1,sprintf('%d active links | %.2f%% margin retained | total weight = %.3f', ...
    nnz(active),100*optimized.margin/rho_ref,sum(w_optimized)), ...
    'Color',gray,'FontSize',12);
exportgraphics(fig1,'Experiment12_Optimized_Topology.png', ...
    'Resolution',300,'BackgroundColor','white');

% Figure 2: fixed-unit versus optimized-weight performance.
fig2 = figure('Color','w','Position',[120 120 1200 650]);
tl = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

ax2 = nexttile(tl);
b2 = bar(ax2,[comparison.Margin_Retained_Percent, ...
              100*comparison.Total_Weight/model.m]);
b2(1).FaceColor = blue;
b2(2).FaceColor = orange;
design_labels = {'Unit weights','Optimized weights'};
set(ax2,'XTick',1:2,'XTickLabel',design_labels, ...
    'XTickLabelRotation',0);
ylabel(ax2,'Percent of reference (%)');
title(ax2,'Stability Retention and Communication Weight','Color','k');
legend(ax2,{'Stability margin retained','Total weight / 6'}, ...
    'Location','northoutside','Orientation','horizontal');
ylim(ax2,[0 110]);
style_publication_axes(ax2);

ax3 = nexttile(tl);
b3 = bar(ax3,[comparison.Common_Settling_Time_s, ...
              comparison.DER4_Settling_Time_s]);
b3(1).FaceColor = blue;
b3(2).FaceColor = orange;
hold(ax3,'on');
yline(ax3,T_common_max,'--','Common limit','Color',blue, ...
    'LineWidth',1.4,'LabelHorizontalAlignment','left');
yline(ax3,T_der4_max,'--','DER 4 limit','Color',orange, ...
    'LineWidth',1.4,'LabelHorizontalAlignment','right');
set(ax3,'XTick',1:2,'XTickLabel',design_labels, ...
    'XTickLabelRotation',0);
ylabel(ax3,'Maximum-node settling time (s)');
title(ax3,'Validated Dynamic Performance','Color','k');
legend(ax3,{'Common input','Localized input at DER 4'}, ...
    'Location','northoutside','Orientation','horizontal');
ylim(ax3,[0 0.50]);
style_publication_axes(ax3);

title(tl,'Experiment 12: Joint Link-and-Weight Optimization', ...
    'Color','k','FontWeight','bold','FontSize',16);
exportgraphics(fig2,'Experiment12_Optimization_Comparison.png', ...
    'Resolution',300,'BackgroundColor','white');

disp('Files saved:');
disp('  Experiment12_Optimization_Starts.csv');
disp('  Experiment12_Optimized_Weights.csv');
disp('  Experiment12_Design_Comparison.csv');
disp('  Experiment12_Optimized_Topology.png');
disp('  Experiment12_Optimization_Comparison.png');

%% Local functions
function J = sparse_objective(w,tau,epsilon_weight)
    smooth_link_count = sum(1-exp(-w/tau));
    J = smooth_link_count + epsilon_weight*sum(w);
end

function [c,ceq] = performance_constraints(w,model,requirements,t_grid)
    metrics = evaluate_design(w,model,t_grid);
    c = [requirements.rho_required-metrics.margin; ...
         metrics.common_settling_time-requirements.T_common_max; ...
         metrics.der4_settling_time-requirements.T_der4_max; ...
         requirements.connectivity_min-metrics.algebraic_connectivity];
    c(~isfinite(c)) = 1e3;
    ceq = [];
end

function metrics = evaluate_design(w,model,t_grid)
    L = weighted_laplacian(model.n,model.possible_edges,w);
    A = -model.h*model.Bel*model.Dp_scaled*L*model.Dp_scaled;
    metrics.margin = effective_stability_margin(A,model.n);

    lap_eigs = sort(real(eig((L+L')/2)));
    metrics.algebraic_connectivity = max(0,lap_eigs(2));

    metrics.common_settling_time = impulse_settling_time(A, ...
        model.B_common,model.capacity,t_grid,model.common_tolerance);
    metrics.der4_settling_time = impulse_settling_time(A, ...
        model.B_der4,model.capacity,t_grid,model.der4_tolerance);
end

function L = weighted_laplacian(n,edges,w)
    W = zeros(n);
    for q = 1:size(edges,1)
        i = edges(q,1); j = edges(q,2);
        W(i,j) = w(q); W(j,i) = w(q);
    end
    L = diag(sum(W,2))-W;
end

function margin = effective_stability_margin(A,n)
    poles = eig(A);
    effective = poles(abs(poles)>1e-7);
    if numel(effective) < n-1 || any(real(effective)>=0)
        margin = 0;
    else
        margin = -max(real(effective));
    end
end

function settling_time = impulse_settling_time(A,B,capacity,t,tolerance)
    % A unit-area impulse initializes the post-disturbance state at B.
    % The response is evaluated directly from the modal solution.
    [V,D] = eig(A);
    lambda = diag(D);
    coefficients = V\B;
    exponentials = exp(lambda*t');
    X = real((V*(coefficients.*exponentials))');
    Xn = X./capacity;
    average_state = mean(Xn,2);
    disagreement = max(abs(Xn-average_state),[],2);

    last_outside = find(disagreement>tolerance,1,'last');
    if isempty(last_outside)
        settling_time = 0;
    elseif last_outside == numel(t)
        settling_time = Inf;
    else
        settling_time = t(last_outside+1);
    end
end

function style_publication_axes(ax)
    ax.Color = 'w';
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.GridColor = [0.75 0.75 0.75];
    ax.GridAlpha = 0.35; ax.FontSize = 11;
    ax.LineWidth = 0.8;
    grid(ax,'on'); box(ax,'on');
    if ~isempty(ax.Toolbar), ax.Toolbar.Visible = 'off'; end
end
