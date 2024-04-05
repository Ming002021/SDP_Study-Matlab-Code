%n is the numbe of nodes
n = 6;
%m is the number of edges
m = 6;
% Define E be the edges of a graph with 6 nodes. Each row in E represents an edge between two nodes.
E = [1 2; 1 3; 2 4; 3 5; 4 6;5 6];
% Define an adjacency matrix A for the graph using the edge information stored in E
A = full(sparse(E(:,1), E(:,2), ones(m,1), n, n)); %[A]_ij=1 means there is an edge between i-j
A = A + A';
s=1;
t=6;
% Test functions

Single_ShortestPath_directedIP(A,s,t);
Single_ShortestPath_directedLP(A,s,t);

% We want to write a function solving the single Source-Sink shortest paths
% problem using IP 

% the input is a connected undirected graph A, the source node s, and the sink node t

function [optimal_path, num_arcs_visited]=Single_ShortestPath_directedIP(A,s,t) %issymmetric(A) 

    % Transform A to be a symmetric matrix
    A = A + A';

    % Get the number of nodes
    n = size(A, 1);

    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, s), t);

    % Define the binary two-dimensional IP decision variable x_{ij}^h, i,j=1..n
    x = binvar(n, n,'full');

    % Objective function: Minimize the total number of arcs visited in the
    % path
  
    obj = sum(x(:));

    % Constraints
    constraints = [];

    % Constraint 1: x_{ij} = 0 if no arc (i,j) exists
    constraints = [constraints, x(:,:) <= A];

    % Constraint 2: Exactly one arc leaving source node s
    constraints = [constraints, sum(x(s,:)) - sum(x(:,s)) == 1];

    % Constraint 3: Exactly one arc entering sink node t
    constraints = [constraints, sum(x(t,:)) - sum(x(:,t)) == -1];

    % Constraint 4: The number of arcs of the path entering v which is not s, t equals the number of arcs of the path leaving v

    for i= N_st
        
        constraints = [constraints,sum(x(i, :)) - sum(x(:, i)) == 0];
    end

    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);

    % Retrieve the optimal solution
    optimal_arc_variables = value(x);

    % Calculate the total number of arcs visited
    num_arcs_visited = sum(optimal_arc_variables(:));

    % Extract the optimal path from the solution
    optimal_path = [];
    current_node = s;
    while current_node ~= t
        next_node = find(optimal_arc_variables(current_node, :));
        optimal_path = [optimal_path, current_node];
        current_node = next_node;
    end
    optimal_path = [optimal_path, t];

end



% We want to write a function solving the single Source-Sink shortest paths
% problem using LP relaxation 

% the input is a connected undirected graph A, the source node s, and the sink node t

function [optimal_path, num_arcs_visited]=Single_ShortestPath_directedLP(A,s,t) %issymmetric(A) 

    % Transform A to be a symmetric matrix
    A = A + A';

    % Get the number of nodes
    n = size(A, 1);

    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, s), t);

    % Define the LP two-dimensional IP decision variable x_{ij}^h, i,j=1..n
    x = sdpvar(n, n,'full');

    % Objective function: Minimize the total number of arcs visited in the
    % path
  
    obj = sum(x(:));

    % Constraints
    constraints = [];
    
   % Constraint 0: x is in [0,1]
    constraints = [constraints, 0 <= x <= 1];

    % Constraint 1: x_{ij} = 0 if no arc (i,j) exists
    constraints = [constraints, x(:,:) <= A];

    % Constraint 2: Exactly one arc leaving source node s
    constraints = [constraints, sum(x(s,:)) - sum(x(:,s)) == 1];

    % Constraint 3: Exactly one arc entering sink node t
    constraints = [constraints, sum(x(t,:)) - sum(x(:,t)) == -1];

    % Constraint 4: The number of arcs of the path entering v not s_h, t_h equals the number of arcs of the path leaving v

    for i= N_st
        
        constraints = [constraints,sum(x(i, :)) - sum(x(:, i)) == 0];
    end

    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);

    % Retrieve the optimal solution
    optimal_arc_variables = value(x);

    % Calculate the total number of arcs visited
    num_arcs_visited = sum(optimal_arc_variables(:));

    % Extract the optimal path from the solution
    optimal_path = [];
    current_node = s;
    while current_node ~= t
        next_node = find(optimal_arc_variables(current_node, :));
        optimal_path = [optimal_path, current_node];
        current_node = next_node;
    end
    optimal_path = [optimal_path, t];

end













