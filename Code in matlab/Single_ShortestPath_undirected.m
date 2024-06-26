%n is the numbe of nodes
n = 6;
%m is the number of edges
m = 6;
% Define E be the edges of a graph with 6 nodes. Each row in E represents an edge between two nodes.
E = [1 2; 1 3; 2 4; 3 5; 4 6;5 6];

% Define an adjacency matrix A for the graph using the edge information stored in E
% Note that A is an upper triangle matrix and we assume for every edge {i,j} in E, we have i < j
A = full(sparse(E(:,1), E(:,2), ones(m,1), n, n)); %[A]_ij=1 means there is an edge between i-j
s=1;
t=6;
% Test functions
Single_ShortestPath_undirectedIP(A,s,t);
Single_ShortestPath_undirectedLP(A,s,t);

% We want to write a function solving the single Source-Sink shortest paths
% problem using IP

% the input is a connected undirected graph A (an upper triangle matrix), the source node s, and the sink node t

function [optimal_path, num_edges_visited]=Single_ShortestPath_undirectedIP(A,s,t) %issymmetric(A) 


    % Get the number of nodes
    n = size(A, 1);

    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, s), t);

    % Define the binary two-dimensional IP decision variable x_{ij}, i,j=1,...n
    x = binvar(n, n,'full');

    % Define the binary one-dimensional IP decision variable y_j,j=1,2,...n
    y = binvar(n, 1,'full');

    % Objective function: Minimize the total number of edges visited in the
    % path
  
    obj = sum(x(:));

    % Constraints
    constraints = [];

    % Constraint 1: x_{ij} = 0 if no edge (i,j) exists. Note the we assume for every edge {i,j} in E, we have i < j
    constraints = [constraints, x(:,:) <= A];

    % Constraint 2: There is exactly one edge incident with s in the path

    constraints = [constraints, sum(x(s,:)) + sum(x(:,s)) == 1];

    % Constraint 3: There is exactly one edge incident with t in the path
    constraints = [constraints, sum(x(t,:)) + sum(x(:,t)) == 1];

    % Constraint 4: The number of edges incident with v which is not s,t, equals 2
    % if this node v is visited

    for i= N_st
        
        constraints = [constraints,sum(x(i, :)) + sum(x(:, i)) == y(i)*2];
    end

    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);

    % Retrieve the optimal solution
    optimal_edge_variables = value(x);
    % Calculate the total number of edges visited
    num_edges_visited = sum(value(obj));

    % Extract the optimal path from the solution
    optimal_path = [];
    current_node = s;
    while current_node ~= t
        next_node = find(optimal_edge_variables(current_node, :));
        optimal_path = [optimal_path, current_node];
        current_node = next_node;
    end
    optimal_path = [optimal_path, t];

    
end


% We want to write a function solving the single Source-Sink shortest paths
% problem using LP relaxation

% the input is a connected undirected graph A (an upper triangle matrix), the source node s, and the sink node t

function [optimal_path, num_edges_visited]=Single_ShortestPath_undirectedLP(A,s,t) %issymmetric(A) 


    % Get the number of nodes
    n = size(A, 1);

    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, s), t);

    % Define the LP two-dimensional IP decision variable x_{ij}, i,j=1,...n
    x = sdpvar(n, n,'full');

    % Define the binary one-dimensional IP decision variable y_j,j=1,2,...n
    y = binvar(n, 1,'full');

    % Objective function: Minimize the total number of edges visited in the
    % path
  
    obj = sum(x(:));

    % Constraints
    constraints = [];
    % Constraint 0: x is in [0,1]
    constraints = [constraints, 0 <= x <= 1];

    % Constraint 1: x_{ij} = 0 if no edge (i,j) exists. Note the we assume for every edge {i,j} in E, we have i < j
    constraints = [constraints, x(:,:) <= A];

    % Constraint 2: There is exactly one edge incident with s in the path

    constraints = [constraints, sum(x(s,:)) + sum(x(:,s)) == 1];

    % Constraint 3: There is exactly one edge incident with t in the path
    constraints = [constraints, sum(x(t,:)) + sum(x(:,t)) == 1];

    % Constraint 4: The number of edges incident with v not s,t, equals 2
    % if this node v is visited

    for i= N_st
        
        constraints = [constraints,sum(x(i, :)) + sum(x(:, i)) == y(i)*2];
    end

    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);

    % Retrieve the optimal solution
    optimal_edge_variables = value(x);
    % Calculate the total number of edges visited
    num_edges_visited = sum(value(obj));

    % Extract the optimal path from the solution
    optimal_path = [];
    current_node = s;
    while current_node ~= t
        next_node = find(optimal_edge_variables(current_node, :));
        optimal_path = [optimal_path, current_node];
        current_node = next_node;
    end
    optimal_path = [optimal_path, t];

    
end











