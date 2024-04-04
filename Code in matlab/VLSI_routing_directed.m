% generate a random graph
n = 22;
p = 1/2;
A = rand(n) > p;
k=5;
s=[13 4 6 2 8];
t=[3 7 15 5 9];

VLSI_routing_directed_IP(A,k,s,t) 


% We want to write a function solving Pairwise Source-Sink vertex-disjoint paths problem

% the input is a connected undirected graph A, and the number of paths we
% want to construct k, the source nodes s, and the sink nodes t

function result = VLSI_routing_directed_IP(A,k,s,t) 
    % Transform A to be a symmetric matrix
    A = A + A';

    % Get the number of nodes
    n = size(A, 1);

    % Get the nodes that are not source nodes and sink nodes
    N_st = setdiff(setdiff(1:n, s), t);

    % Define the binary three-dimensional IP decision variable x_{ij}^h, i,j=1..n
    x = binvar(n, n, k, 'full');

    % Objective function: Minimize the total number of arcs visited in all paths
    obj = sum(x(:));

    % Constraints
    constraints = [];
    for h = 1:k
        % Constraint 1: x_{ij}^h = 0 if no arc (i,j) exists
        constraints = [constraints, x(:,:,h) <= A];

        % Constraint 2: Exactly one arc leaving source node si
        constraints = [constraints, sum(x(s(h), :, h)) - sum(x(:, s(h), h)) == 1];

        % Constraint 3: Exactly one arc entering sink node ti
        constraints = [constraints, sum(x(t(h), :, h)) - sum(x(:, t(h), h)) == -1];

        % Constraint 4: In each path Ph, the number of arcs of the path entering v not s_h, t_h equals the number of arcs of the path leaving v

        for i= N_st
            
            constraints = [constraints,sum(x(i, :, h)) - sum(x(:, i, h)) == 0];
        end
        
    end

    % Vertex disjoint constraints
    for i = N_st
        constraints = [constraints, sum(sum(x(i, :, :))) <= 1, sum(sum(x(:, i, :))) <= 1];
    end

    % Solve the optimization problem
    options = sdpsettings('solver', 'gurobi'); % Use the solver gurobi
    optimize(constraints, obj, options);

    disp(value(x));

    % Retrieve the optimal objective value
    result = value(obj);

end
