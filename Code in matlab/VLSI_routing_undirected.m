% generate a random graph
n = 22;
p = 1/2;
A = rand(n) > p;
% Note that A should be an upper triangle matrix and we assume for every edge {i,j} in E, we have i < j
A = triu(A,1);
s=[13 4 6 2 8];
t=[3 7 15 5 9];
VLSI_routing_undirected_IP(A,s,t);
VLSI_routing_undirected_LP(A,s,t)

% We want to write a function solving Pairwise Source-Sink vertex-disjoint
% paths problem using IP

% the input is a connected undirected graph A, the source nodes s, and the sink nodes t

function [optimal_paths, num_arcs_visited] = VLSI_routing_undirected_IP(A, s, t) 
    
    % Get the number of nodes
    n = size(A, 1);
    
    % Get the number of paths
    k = size(s, 2);
    
    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, unique(s)), unique(t));
    
    % Define the binary three-dimensional IP decision variable x_{ij}^h, i,j=1..n
    x = binvar(n, n, k, 'full');
    
    % Define the binary two-dimensional IP decision variable y_j^h, j=1,2,...,n
    y = binvar(n, k, 'full');
    
    % Define the objective function: Minimize the total number of arcs visited in all paths
    obj = sum(x(:));
    
    % Constraints
    constraints = [];
    
    for h = 1:k
        % Constraint 1: x_{ij}^h = 0 if no edge (i,j) exists. Note the we assume for every edge {i,j} in E, we have i < j
        constraints = [constraints, x(:,:,h) <= A];
    
        % Constraint 2: Exactly one edge incident with source node s(h)
        constraints = [constraints, sum(x(s(h), :, h)) + sum(x(:, s(h), h)) == 1];
    
        % Constraint 3: Exactly one edge incident with sink node t(h)
        constraints = [constraints, sum(x(t(h), :, h)) + sum(x(:, t(h), h)) == 1];
    
        % Constraint 4: The number of edges incident with node v, not s(h), t(h), equals 2
        for i = N_st
            constraints = [constraints, sum(x(i, :, h)) + sum(x(:, i, h)) == y(i,h)*2];
        end
    end
    
    % Vertex disjoint constraints
    for i = N_st
        constraints = [constraints, sum(sum(x(i, :, :))) <= 1, sum(sum(x(:, i, :))) <= 1];
    end
    
    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);
    
    % Retrieve the optimal solution
    optimal_arc_variables = value(x);
    
    % Calculate the total number of arcs visited
    num_arcs_visited = sum(optimal_arc_variables(:));
    
    % Extract the optimal paths from the solution
    optimal_paths = cell(1, k);
    
    for h = 1:k
        optimal_path_h = [];
        current_node = s(h);
        while current_node ~= t(h)
            next_node = find(optimal_arc_variables(current_node, :, h));
            optimal_path_h = [optimal_path_h, current_node];
            current_node = next_node;
        end
        optimal_path_h = [optimal_path_h, t(h)];
        optimal_paths{h} = optimal_path_h;
    end
end


% We want to write a function solving Pairwise Source-Sink vertex-disjoint
% paths problem using LP relaxtion

% the input is a connected undirected graph A, the source nodes s, and the sink nodes t

function [optimal_paths, num_arcs_visited] = VLSI_routing_undirected_LP(A, s, t) 
    
    % Get the number of nodes
    n = size(A, 1);
    
    % Get the number of paths
    k = size(s, 2);
    
    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, unique(s)), unique(t));
    
    % Define the LP three-dimensional IP decision variable x_{ij}^h, i,j=1..n
    x = sdpvar(n, n, k, 'full');
    
    % Define the binary two-dimensional IP decision variable y_j^h, j=1,2,...,n
    y = binvar(n, k, 'full');
    
    % Define the objective function: Minimize the total number of arcs visited in all paths
    obj = sum(x(:));
    
    % Constraints
    constraints = [];
    
    % Constraint 0: x is in [0,1]
    constraints = [constraints, 0 <= x <= 1];

    for h = 1:k
        % Constraint 1: x_{ij}^h = 0 if no edge (i,j) exists. Note the we assume for every edge {i,j} in E, we have i < j
        constraints = [constraints, x(:,:,h) <= A];
    
        % Constraint 2: Exactly one edge incident with source node s(h)
        constraints = [constraints, sum(x(s(h), :, h)) + sum(x(:, s(h), h)) == 1];
    
        % Constraint 3: Exactly one edge incident with sink node t(h)
        constraints = [constraints, sum(x(t(h), :, h)) + sum(x(:, t(h), h)) == 1];
    
        % Constraint 4: The number of edges incident with node v, not s(h), t(h), equals 2
        for i = N_st
            constraints = [constraints, sum(x(i, :, h)) + sum(x(:, i, h)) == y(i,h)*2];
        end
    end
    
    % Vertex disjoint constraints
    for i = N_st
        constraints = [constraints, sum(sum(x(i, :, :))) <= 1, sum(sum(x(:, i, :))) <= 1];
    end
    
    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);
    
    % Retrieve the optimal solution
    optimal_arc_variables = value(x);
    
    % Calculate the total number of arcs visited
    num_arcs_visited = sum(optimal_arc_variables(:));
    
    % Extract the optimal paths from the solution
    optimal_paths = cell(1, k);
    
    for h = 1:k
        optimal_path_h = [];
        current_node = s(h);
        while current_node ~= t(h)
            next_node = find(optimal_arc_variables(current_node, :, h));
            optimal_path_h = [optimal_path_h, current_node];
            current_node = next_node;
        end
        optimal_path_h = [optimal_path_h, t(h)];
        optimal_paths{h} = optimal_path_h;
    end
end

