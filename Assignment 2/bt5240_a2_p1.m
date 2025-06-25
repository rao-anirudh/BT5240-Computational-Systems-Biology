% BT5240 - Assignment 2 - Problem 1
% Anirudh Rao (BE21B004)

%% Preliminaries

% Adding MATLAB BGL to the PATH

addpath(genpath('matlab_bgl-4.0.1'))

% Clearing screen and variables

clc, clearvars;

% Setting a random seed for replicability

rng(5240)

% NOTE: Ensure that MATLAB's Statistics and Machine Learning Toolbox has been installed

%% Problem 1

% Function to create a Watts-Strogatz network with n nodes, k initial neighbours, and p rewiring probability

function A = watts_strogatz(n, k, p)

    % Initialising the adjacency matrix of the graph

    A = zeros(n);

    % Creating a ring lattice where each of the n nodes is connected to k of its neighbours

    for neighbour = 1:k/2
        A = A + circshift(eye(n), neighbour) + circshift(eye(n), -neighbour);
    end

    % Rewiring the edges with probability p

    A_new = A; 
    [i, j] = find(triu(A)); 
    edges = [i j]; 
    
    for e = 1:size(edges, 1)

        node = edges(e, 1);
        neighbour = edges(e, 2);
        
        if rand < p
            
            A_new(node, neighbour) = 0;
            A_new(neighbour, node) = 0;
           
            possible_neighbours = setdiff(1:n, [node; find(A_new(node, :))']); 
            
            if ~isempty(possible_neighbours) 
                new_neighbour = possible_neighbours(randi(length(possible_neighbours)));
                A_new(node, new_neighbour) = 1;
                A_new(new_neighbour, node) = 1;
            else
                A_new(node, neighbour) = 1;
                A_new(neighbour, node) = 1;
            end
            
        end

    end

    A = sparse(A_new);

end

% Parameters given in the question

n = 100;
k = 10;
P = [0 0.3 0.7 1];
edge_prob = ((n*k)/2) / (n*(n-1)/2);

for p = P
    
    % Constructing different Watts-Strogatz networks

    A = watts_strogatz(n, k, p);
    G = graph(A);
    figure()
    plot(G, 'Layout', 'circle')
    if p == 0
        title(sprintf('Ring Lattice (n=%d, k=%d)', n, k));
    else
        title(sprintf('Watts-Strogatz Network (n=%d, k=%d, p=%.2f)', n, k, p));
    end
    
    if p ~= 0
        fprintf('<strong>Watts-Strogatz network with (n=%d, k=%d, p=%.2f)</strong>\n', n, k, p)
    
        % Computing and plotting the degree distribution
    
        deg = degree(G);
        figure()
        histogram(deg)
        xlabel('Degree');
        ylabel('Frequency');
        title(sprintf('Watts-Strogatz Degree Distribution (n=%d, k=%d, p=%.2f)', n, k, p));
        
        % Computing the average clustering coefficient
    
        average_clustering = mean(clustering_coefficients(A));
    
        % Computing the characteristic path length
        
        dists = all_shortest_paths(A);
        dists = triu(dists);
        char_path_length = mean(dists,'all');
        
        % Computing the properties of 100 random networks
    
        num_random = 100;
    
        random_average_clusterings = zeros(num_random, 1);
        random_char_path_lengths = zeros(num_random, 1);
    
        for i = 1:num_random
    
            % Constructing an Erdos-Renyi random network
    
            A_random = erdos_reyni(n, edge_prob);
    
            % Computing the average clustering coefficient of the random network
    
            random_average_clusterings(i) = mean(clustering_coefficients(A_random));
    
            % Computing the characteristic path length of the random network
            
            random_dists = all_shortest_paths(A_random);
            random_dists = triu(random_dists);
            if mean(random_dists, 'all') ~= Inf
                random_char_path_lengths(i) = mean(random_dists,'all');
            end
    
        end
    
        random_char_path_lengths = random_char_path_lengths(random_char_path_lengths ~= 0);
        
        % Performing statistical tests
      
        random_mean_clustering = mean(random_average_clusterings);
        random_std_clustering = std(random_average_clusterings);
        random_mean_path_length = mean(random_char_path_lengths);
        random_std_path_length = std(random_char_path_lengths);
        
        [~,pval,~,zval] = ztest(average_clustering, random_mean_clustering, random_std_clustering);
        fprintf('Average clustering coefficient = %f \n', average_clustering)
        fprintf('Z-score for average clustering coefficient = %f\n', zval)
        fprintf('p-value for average clustering coefficient = %f\n', pval)

        figure()
        histfit(random_average_clusterings);
        hold on;
        y_limits = ylim;
        plot([average_clustering, average_clustering], y_limits, 'r', 'LineWidth', 2);
        xlabel('Average Clustering Coefficient');
        ylabel('Frequency');
        title(sprintf('Watts-Strogatz (n=%d, k=%d, p=%.2f) - Average Clustering Coefficient', n, k, p));
        legend('', 'Erdős–Rényi', 'Watts-Strogatz', 'Location', 'best');
        hold off;
        
        [~,pval,~,zval] = ztest(char_path_length, random_mean_path_length, random_std_path_length);
        fprintf('Characteristic path length = %f \n', char_path_length)
        fprintf('Z-score for characteristic path length = %f\n', zval)
        fprintf('p-value for characteristic path length = %f\n\n', pval)

        figure()
        histfit(random_char_path_lengths);
        hold on;
        y_limits = ylim;
        plot([char_path_length, char_path_length], y_limits, 'r', 'LineWidth', 2);
        xlabel('Characteristic Path Length');
        ylabel('Frequency');
        title(sprintf('Watts-Strogatz (n=%d, k=%d, p=%.2f) - Characteristic Path Length', n, k, p));
        legend('', 'Erdős–Rényi', 'Watts-Strogatz', 'Location', 'best');
        hold off;

    end

end
