%% EE 800 Experiment 6: four-node communication-topology enumeration
% Enumerates every connected, unweighted, undirected graph on four labeled
% DER nodes. The physical model, capacities, gain, and disturbance are held
% fixed. O-NAPC and A-NAPC are compared by effective stability margin.

clc;
clear;
close all;

%% Fixed model from Experiments 1 through 5
n = 4;
h = 0.1;

Bel = [ 10 -10   0   0;
       -10  20 -10   0;
         0 -10  20 -10;
         0   0 -10  10];

capacity = [8 4 2 1];
Dp_scaled = diag([1 2 4 8]);

% All six possible undirected edges for four labeled nodes.
edges = [1 2;
         1 3;
         1 4;
         2 3;
         2 4;
         3 4];

%% Enumerate all nonempty edge subsets and retain connected graphs
records = struct('Graph_ID',{},'Mask',{},'Link_Count',{},'Edges',{}, ...
    'O_NAPC_Margin',{},'A_NAPC_Margin',{}, ...
    'O_Dominant_Real',{},'O_Dominant_Imag',{}, ...
    'A_Dominant_Real',{},'A_Dominant_Imag',{});

graph_id = 0;

for mask = 1:(2^size(edges,1)-1)
    active = logical(bitget(mask,1:size(edges,1)));
    selected_edges = edges(active,:);

    % Build the communication adjacency matrix and Laplacian.
    W = zeros(n);
    for e = 1:size(selected_edges,1)
        i = selected_edges(e,1);
        j = selected_edges(e,2);
        W(i,j) = 1;
        W(j,i) = 1;
    end
    Lcom = diag(sum(W,2)) - W;

    % For an undirected four-node graph, rank(L) = n-1 means connected.
    if rank(Lcom,1e-10) ~= n-1
        continue;
    end

    graph_id = graph_id + 1;

    % Same controller matrices used in the earlier controlled experiments.
    A_O = -h * Bel * Lcom * Dp_scaled;
    A_A = -h * Bel * Dp_scaled * Lcom * Dp_scaled;

    [margin_O,dominant_O] = effectiveMargin(A_O);
    [margin_A,dominant_A] = effectiveMargin(A_A);

    edge_text = strings(size(selected_edges,1),1);
    for e = 1:size(selected_edges,1)
        edge_text(e) = sprintf('%d-%d',selected_edges(e,1),selected_edges(e,2));
    end

    records(graph_id).Graph_ID = graph_id;
    records(graph_id).Mask = mask;
    records(graph_id).Link_Count = size(selected_edges,1);
    records(graph_id).Edges = strjoin(edge_text,', ');
    records(graph_id).O_NAPC_Margin = margin_O;
    records(graph_id).A_NAPC_Margin = margin_A;
    records(graph_id).O_Dominant_Real = real(dominant_O);
    records(graph_id).O_Dominant_Imag = imag(dominant_O);
    records(graph_id).A_Dominant_Real = real(dominant_A);
    records(graph_id).A_Dominant_Imag = imag(dominant_A);
end

Results = struct2table(records);
Results = sortrows(Results,{'Link_Count','A_NAPC_Margin'}, ...
    {'ascend','descend'});

fprintf('Connected labeled graphs found: %d\n',height(Results));
disp(Results);
writetable(Results,'Experiment6_All_Connected_Topologies.csv');

%% Best graph for each feasible link count
link_counts = unique(Results.Link_Count);
best_rows = zeros(numel(link_counts),1);

for k = 1:numel(link_counts)
    candidates = find(Results.Link_Count == link_counts(k));
    [~,local_index] = max(Results.A_NAPC_Margin(candidates));
    best_rows(k) = candidates(local_index);
end

BestByLinkCount = Results(best_rows,:);
fprintf('\nBEST A-NAPC TOPOLOGY AT EACH LINK COUNT\n');
disp(BestByLinkCount(:,{'Link_Count','Graph_ID','Edges', ...
    'O_NAPC_Margin','A_NAPC_Margin'}));
writetable(BestByLinkCount,'Experiment6_Best_By_Link_Count.csv');

%% Figure 1: every connected topology
figure('Name','All connected communication topologies');
scatter(Results.Link_Count,Results.O_NAPC_Margin,55,'filled');
hold on;
scatter(Results.Link_Count,Results.A_NAPC_Margin,55,'filled');
grid on;
xlabel('Number of undirected communication links');
ylabel('Effective stability margin');
title('All Connected Four-Node Communication Topologies');
legend({'O-NAPC','A-NAPC'},'Location','northwest');
xticks(3:6);
exportgraphics(gcf,'Experiment6_All_Topologies.png', ...
    'Resolution',300,'BackgroundColor','current');

%% Figure 2: best A-NAPC design at each communication cost
figure('Name','Best margin by link count');
plot(BestByLinkCount.Link_Count,BestByLinkCount.O_NAPC_Margin, ...
    '-o','LineWidth',1.6,'MarkerSize',7);
hold on;
plot(BestByLinkCount.Link_Count,BestByLinkCount.A_NAPC_Margin, ...
    '-s','LineWidth',1.6,'MarkerSize',7);
grid on;
xlabel('Number of undirected communication links');
ylabel('Effective stability margin');
title('Best A-NAPC Topology at Each Link Count');
legend({'O-NAPC margin of selected graph','A-NAPC margin'}, ...
    'Location','northwest');
xticks(link_counts);
exportgraphics(gcf,'Experiment6_Best_By_Link_Count.png', ...
    'Resolution',300,'BackgroundColor','current');

%% Figure 3: draw the selected best communication graphs
figure('Name','Selected communication graphs');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Fixed node positions:
% DER 1 = right, DER 2 = top, DER 3 = left, DER 4 = bottom
x_position = [1 0 -1 0];
y_position = [0 1 0 -1];

for k = 1:height(BestByLinkCount)

    mask = BestByLinkCount.Mask(k);
    active = logical(bitget(mask,1:size(edges,1)));
    selected_edges = edges(active,:);

    G = graph(selected_edges(:,1),selected_edges(:,2),[],n);

    nexttile;

    % Plot nodes and edges without MATLAB's automatic labels
    topology_plot = plot(G, ...
        'XData',x_position, ...
        'YData',y_position, ...
        'NodeLabel',{}, ...
        'LineWidth',1.5, ...
        'MarkerSize',8);

    hold on;

    % Add labels manually so they remain horizontal
    % Labels positioned above or below each node
    text(1,0.14,'DER 1', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'Color','white');

    text(0,1.14,'DER 2', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'Color','white');

    text(-1,0.14,'DER 3', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'Color','white');

    text(0,-1.14,'DER 4', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'Color','white');

    title(sprintf('%d links, A-margin = %.4f', ...
        BestByLinkCount.Link_Count(k), ...
        BestByLinkCount.A_NAPC_Margin(k)));

    xlim([-1.40 1.40]);
    ylim([-1.35 1.35]);
    pbaspect([1 1 1]);

    set(gca, ...
        'XTick',[], ...
        'YTick',[], ...
        'Box','on', ...
        'XColor',[0.45 0.45 0.45], ...
        'YColor',[0.45 0.45 0.45], ...
        'LineWidth',0.8);
end

exportgraphics(gcf,'Experiment6_Selected_Topologies.png', ...
    'Resolution',300);

%% Local function: remove the expected consensus eigenvalue
function [margin,dominant_pole] = effectiveMargin(A)
    lambda = eig(A);
    [~,zero_index] = min(abs(lambda));
    lambda_effective = lambda;
    lambda_effective(zero_index) = [];

    [dominant_real,index] = max(real(lambda_effective));
    dominant_pole = lambda_effective(index);
    margin = -dominant_real;
end
