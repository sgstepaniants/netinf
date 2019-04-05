function Y = harmonicOsc(K, p, v, m, c, tSpan, bc, forcingFunc)
% calls ode45 to solve nearest neighbor coupled spring mass model
% 
% INPUTS:
%
% K
%       a matrix of spring coefficents (the underlying connectivity matrix)
%
% p
%       an [n x 1] vector of initial positions for the ODE
%
% v
%       an [n x 1] vector of initial velocities for the ODE
%
% adj
%       an [n x n] adjacency matrix for the network
%
% m 
%       an [n x 1] vector of masses, one for each node
%       in the network
%
% c 
%       an n - 1, n, or n + 1 length vector of damping coefficients,
%       one for each edge in the network
%
% tSpan
%       a [1 x m] vector of times for which you want the solution of the
%       ODE returned
%
% bc
%       a string that specifies the boundary conditions of the nearest
%       neighbor spring model ('circ', 'fixed', 'free')
%
% pert  a number that signifies how much nodes should be perturbed when
%       their oscillation amplitudes fall below a certain threshold
%
% OUTPUTS:
%
% Y
%       an [n x m] matrix of the solution of the ODE for all n 
%       nodes and m time points

    n = size(p, 1);
    param{1} = n;     % number of nodes
    param{2} = K;     % spring constants
    param{3} = m;     % masses
    param{4} = c;     % damping coefficients
    param{5} = bc;    % boundary conditions
    
    if nargin < 8
        forcingFunc = zeros(n, length(tSpan));
    end
    
    param{6} = @(t) forcingFunc(:, findClosestTime(t, tSpan));
    
    y0 = [p; v];
    
    % Solve this ode
    SOL = ode15s(@(t, y) odeHarmonic(t, y, param), tSpan, y0);
    Y = deval(SOL, tSpan);
    
    % Take only the top portion of the matrix which holds the positions
    % of each one of the n masses over time
    Y = Y(1:n, :, :);
end

function idx = findClosestTime(t, tSpan)
    [~, idx] = min(abs(tSpan-t));
end