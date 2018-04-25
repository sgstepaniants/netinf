function Y = GenerateNNCoupledData(n, tSpan, N, randpfn, randvfn, randmfn, randkfn, randcfn, bc, pert)
%
% calls nncoupled.m to return solution of NNCoupled model
%
% INPUTS:
%
% n
%       number of nodes in the network
%
% tSpan
%       a [1 x m] vector of times for which you want the solution of the
%       ODE returned
%
% N
%       a scalar value for the number of random trials requested 
%
%
% randpfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial positions, one for each node in the 
%       network
%
% randvfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial velocities, one for each node in the 
%       network
%
% randmfn 
%       a function that takes a scalar n and retuns an [n x 1] vector 
%       of (possibly random) masses, one for each node
%       in the network
%
% randkfn 
%       a function that takes a scalar n and retuns a vector of length
%       n - 1, n, or n + 1 of (possibly random) spring constants,
%       one for each node in the network
%
% randcfn 
%       a function that takes a scalar n and retuns a vector of length
%       n - 1, n, or n + 1 of (possibly random) damping coefficients,
%       one for each node in the network
%
% bc
%       a string that specifies the boundary conditions of the nearest
%       neighbor spring model ('circ', 'fixed', 'free')
%
% pert  a number that signifies how much nodes should be perturbed when
%       their velocities and accelerations fall below a certain threshold
%
% OUTPUTS:
%
% Y
%       an [n x m x N] matrix of the solution of the ODE for n nodes,
%       m time points, and N random trials

    if nargin == 8
        bc = 'circ';
        pert = 0;
    elseif nargin == 9
        pert = 0;
    end

    m = length(tSpan);
    Y = zeros(n, m, N);
    
    for j = 1 : N
        pos = randpfn(n);
        vel = randvfn(n);
        
        mass = randmfn(n);
        spring = randkfn(n);
        damping = randcfn(n);
        
        if size(spring, 1) < n - 1 || size(spring, 1) > n + 1
            error('spring constants must have length n - 1, n, or n + 1')
        end
        
        if size(damping, 1) < n - 1 || size(damping, 1) > n + 1
            error('damping coefficients must have length n - 1, n, or n + 1')
        end
        
        Y(:, :, j) = nncoupled(pos, vel, mass, spring, damping, tSpan, bc, pert);
    end
    
    % if N == 1, want 2D Y
    Y = squeeze(Y);
end
