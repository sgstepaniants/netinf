% Based on the differential equation for the nncoupled system
% [x', v'] = A[x, v], create A matrix and find its eigenvalues. Then stack
% this system into a Hankel matrix s times and compute its eigenvalues.
% In this demo, we assume all masses are 1.

% number of nodes
n = 10;
% spring constant (same for all nodes)
k = 1;
% damping coefficient (same for all nodes)
c = 1;

E = zeros(n);
F = eye(n);

v = k * ones([1, n-1]);
G = diag(v, 1) - diag(v, -1) - 2 * k * eye(n);
G(1, n) = -k;
G(n, 1) = k;

J = -c * eye(n);

% create the block diagonal matrix A
A = [[E, F]; [G, J]];

% create the Hankel matrix for our differential equation
s = 1;  % number of stacks
