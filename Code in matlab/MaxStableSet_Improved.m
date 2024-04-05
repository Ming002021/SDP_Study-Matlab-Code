%n is the numbe of nodes
n = 6;
%m is the number of edges
m = 10;
% Define E be the edges of a graph with 6 nodes. Each row in E represents an edge between two nodes.
E = [1 2; 1 4; 1 5; 1 6; 2 4; 2 5; 3 4; 3 6; 4 5; 5 6];
% Define an adjacency matrix A for the graph using the edge information stored in E
A = full(sparse(E(:,1), E(:,2), ones(m,1), n, n)); %[A]_ij=1 means there is an edge between i,j
A = A + A';

% generate a random graph
n = 22;
p = 1/2;
A = rand(n) > p;
A = triu(A,1) + triu(A,1)';
%%% 1) We use SDP relaxation to solve it
% Define an SDP variable X of size (n+1) x (n+1)
X = sdpvar(n+1);
% Constraints
% Constraints: x_i * x_j =0 for all adjacent nodes
% tril(blkdiag(0,A),-1) outputs extract only the elements below the main diagonal.
con_1 = [X(find(tril(blkdiag(0,A),-1))) == 0];
con_SDP = [con_1, X >= 0, X(1,1) == 1, X(1, 2:end)' == diag(X(2:end, 2:end))];
obj_SDP = trace(X)-1;
optimize(con_SDP, -obj_SDP)
value(obj_SDP)
%% The output is 'Successfully solved (MOSEK-SDP)'
opt_solution_SDP=value(X);
opt_value_SDP=value(obj_SDP);

%%% 2) We use LP relaxation to solve it
% Define LP decidion variable
x = sdpvar(n, 1);
% Objective function: Maximize the cut
obj = sum(x);
% Constraints: x_i + x_j <= 1 for all adjacent nodes
con = [];
for i = 1:n
for j = i+1:n
if A(i,j) == 1 % If there is an edge between nodes i and j
con = [con; x(i) + x(j) <= 1];
end
end
end
% Constraints: x_i in [0, 1]
con_LP = [con, 0 <= x <= 1];
% Maximize the objective function
optimize(con_LP, -obj)
%% The output is 'Successfully solved (MOSEK-LP/QP)'
opt_solution_LP=value(x);
opt_value_LP=value (obj);
%%% 3) We use IP to solve it
% Constraints: x_i in {0, 1}
con_IP = [con, 0 <= x <= 1, binary(x)];
% Maximize the objective function
optimize(con_IP, -obj)
%% The output is 'Successfully solved (MOSEK-LP/QP-INTEGER)'
opt_solution_IP=value(x);
opt_value_IP=value (obj);
disp(opt_value_IP) % 2
disp(opt_value_LP) % 3
disp(opt_value_SDP) % 2