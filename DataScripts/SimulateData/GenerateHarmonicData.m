function Y = GenerateHarmonicData(n, tSpan, N, K, pfn, vfn, mfn, cfn, bc, forcingFunc)
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
% K
%       a matrix of spring coefficents (the underlying connectivity matrix)
%
% pfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial positions, one for each node in the 
%       network
%
% vfn 
%       a function that takes a scalar n and returns an [n x 1] vector 
%       of (possibly random) initial velocities, one for each node in the 
%       network
%
% mfn 
%       a function that takes a scalar n and retuns an [n x 1] vector 
%       of (possibly random) masses, one for each node
%       in the network
%
%
% cfn 
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

    if nargin == 9
        forcingFunc = zeros(n, length(tSpan));
    end

    m = length(tSpan);
    Y = zeros(n, m, N);
    
    for j = 1 : N
        pos = pfn(n);
        vel = vfn(n);
        
        mass = mfn(n);
        damping = cfn(n);
        
        if strcmp(bc, 'circ')
            if size(K, 1) ~= n
                error('spring constants for circular boundary conditions must have dimension n x n')
            end
        elseif strcmp(bc, 'fixed')
            if size(K, 1) ~= n + 2
                error('spring constants for fixed boundary conditions must have dimension (n + 2) x (n + 2)')
            end
        elseif strcmp(bc, 'free')
            if size(K, 1) ~= n
                error('spring constants for free boundary conditions must have dimension n x n')
            end
        else
            error('invalid boundary conditions')
        end
        
        
        Y(:, :, j) = harmonicOsc(K, pos, vel, mass, damping, tSpan, bc, forcingFunc);
    end
    
    % if N == 1, want 2D Y
    Y = squeeze(Y);
end
