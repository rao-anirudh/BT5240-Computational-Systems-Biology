% BT5240 - Assignment 2 - Problem 2
% Anirudh Rao (BE21B004)

%% Preliminaries

% Adding MATLAB BGL to the PATH

addpath(genpath('matlab_bgl-4.0.1'))

% Clearing screen and variables

clc, clearvars;

%% Problem 2a

% Reading the data

data = readtable('Assignment2_2.txt','FileType','text');

% Extracting the operons, transcription factors, and regulation types

operons = string(data.operon);
tfs = string(data.transcription_factor);
reg_types = string(data.regulation_type);

% Creating a directed graph

G = digraph(tfs, operons);
A = adjacency(G);

% Plotting the network

figure();
graph_plot = plot(G, 'Layout', 'force'); 
title('E. coli Transcription Network');

% Labelling the nodes

labelnode(graph_plot, 1:numnodes(G), G.Nodes.Name);

% Defining edge colours for different regulation types

repressor_colour = [1, 0, 0]; % Red
activator_colour = [0, 1, 0]; % Green
dual_colour = [0, 0, 1]; % Blue

% Looping through the edges and assigning colours individually

for i = 1:height(data)

    source_node = data.transcription_factor(i);
    target_node = data.operon(i);
    
    if data.regulation_type(i) == "activator"
        edge_colour = activator_colour;
    elseif data.regulation_type(i) == "repressor"
        edge_colour = repressor_colour;
    else
        edge_colour = dual_colour;
    end
    
    highlight(graph_plot, source_node, target_node, 'EdgeColor', edge_colour);

end

%% Problem 2b

fprintf('<strong>COMPLETE NETWORK</strong>\n\n');

% Computing the degree centrality and finding the top 5 nodes

deg = centrality(G, 'indegree') + centrality(G, 'outdegree');
[degree, idx] = maxk(deg, 5);
tf_name = G.Nodes.Name(idx);
degree_centrality = int32(degree);
fprintf('<strong>Top 5 TFs by Degree Centrality</strong>\n');
disp(table(tf_name, degree_centrality));

% Computing the closeness centrality and finding the top 5 nodes

closeness = centrality(G, 'incloseness') + centrality(G, 'outcloseness');
[closeness_centrality, idx] = maxk(closeness, 5);
tf_name = G.Nodes.Name(idx);
fprintf('<strong>Top 5 TFs by Closeness Centrality</strong>\n');
disp(table(tf_name, closeness_centrality));

% Computing the betweenness centrality and finding the top 5 nodes

betweenness = centrality(G, 'betweenness');
[betweenness_centrality, idx] = maxk(betweenness, 5);
tf_name = G.Nodes.Name(idx);
fprintf('<strong>Top 5 TFs by Betweenness Centrality</strong>\n');
disp(table(tf_name, betweenness_centrality));

%% Problem 2c

fprintf('<strong>ACTIVATOR NETWORK</strong>\n\n');

% Computing the 'activation fraction' of each of the transcription factors

total_network_activation = sum(strcmp(data.regulation_type, 'activator'));
activation_fraction = zeros(size(G.Nodes, 1), 1);

for i = 1:size(G.Nodes, 1)
    tf = G.Nodes.Name(i);
    tf_activation = sum(strcmp(data.transcription_factor, tf) & strcmp(data.regulation_type, 'activator'));
    activation_fraction(i) = tf_activation / total_network_activation; 
end

% Finding the top 5 activators

[activation_fraction, idx] = maxk(activation_fraction, 5);
tf_name = G.Nodes.Name(idx);
fprintf('<strong>Top 5 Activators</strong>\n');
disp(table(tf_name, activation_fraction));
fprintf('These activators will be removed from the network. Only activating edges will be allowed to remain in the network.\n\n')

% Removing the top 5 activators and retaining only activations in the network

filtered_data = data(~ismember(data.transcription_factor, tf_name), :);
filtered_data = filtered_data(strcmp(filtered_data.regulation_type, 'activator'), :);
G_activator = digraph(filtered_data.transcription_factor, filtered_data.operon);

% Plotting the degree distribution, computing the degree centrality and finding the top 5 nodes

deg = centrality(G_activator, 'indegree') + centrality(G_activator, 'outdegree');

figure()
histogram(deg)
xlabel('Degree');
ylabel('Frequency');
title('Degree Distribution of Activator Network')

[degree, idx] = maxk(deg, 5);
tf_name = G_activator.Nodes.Name(idx);
degree_centrality = int32(degree);
fprintf('<strong>Top 5 Activators by Degree Centrality</strong>\n');
disp(table(tf_name, degree_centrality));

% Computing the closeness centrality and finding the top 5 nodes

closeness = centrality(G_activator, 'incloseness') + centrality(G_activator, 'outcloseness');
[closeness_centrality, idx] = maxk(closeness, 5);
tf_name = G_activator.Nodes.Name(idx);
fprintf('<strong>Top 5 Activators by Closeness Centrality</strong>\n');
disp(table(tf_name, closeness_centrality));

%% Problem 2d

fprintf('<strong>REPRESSOR/DUAL REGULATION NETWORK</strong>\n\n');

% Computing the 'repression fraction' of each of the transcription factors

total_network_repression = sum(strcmp(data.regulation_type, 'repressor'));
repression_fraction = zeros(size(G.Nodes, 1), 1);

for i = 1:size(G.Nodes, 1)
    tf = G.Nodes.Name(i);
    tf_repression = sum(strcmp(data.transcription_factor, tf) & strcmp(data.regulation_type, 'repressor'));
    repression_fraction(i) = tf_repression / total_network_repression; 
end

% Finding the top 5 repressors

[repression_fraction, idx] = maxk(repression_fraction, 5);
tf_name = G.Nodes.Name(idx);
fprintf('<strong>Top 5 Repressors</strong>\n');
disp(table(tf_name, repression_fraction));
fprintf('These repressors will be removed from the network. Only repressing and dual edges will be allowed to remain in the network.\n\n')

% Removing the top 5 repressors and retaining only repressions/dual regulations in the network

filtered_data = data(~ismember(data.transcription_factor, tf_name), :);
filtered_data = filtered_data(strcmp(filtered_data.regulation_type, 'repressor') | strcmp(filtered_data.regulation_type, 'dual'), :);
G_repressor_dual = digraph(filtered_data.transcription_factor, filtered_data.operon);

% Plotting the degree distribution, computing the degree centrality and finding the top 5 nodes

deg = centrality(G_repressor_dual, 'indegree') + centrality(G_repressor_dual, 'outdegree');

figure()
histogram(deg)
xlabel('Degree');
ylabel('Frequency');
title('Degree Distribution of Repressor/Dual Regulation Network')

[degree, idx] = maxk(deg, 5);
tf_name = G_repressor_dual.Nodes.Name(idx);
degree_centrality = int32(degree);
fprintf('<strong>Top 5 Repressors/Dual TFs by Degree Centrality</strong>\n');
disp(table(tf_name, degree_centrality));

% Computing the closeness centrality and finding the top 5 nodes

closeness = centrality(G_repressor_dual, 'incloseness') + centrality(G_repressor_dual, 'outcloseness');
[closeness_centrality, idx] = maxk(closeness, 5);
tf_name = G_repressor_dual.Nodes.Name(idx);
fprintf('<strong>Top 5 Repressors/Dual TFs by Closeness Centrality</strong>\n');
disp(table(tf_name, closeness_centrality));
